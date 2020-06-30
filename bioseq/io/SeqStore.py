import gzip
import os
import sys
import subprocess as sp

import Bio.SeqIO as bpio
from config.settings.base import SEQSTORE_DIR_TMP, SEQSTORE_DIR
from Bio import SeqIO

from SNDG import Struct
import tempfile


class SeqStream:

    def __init__(self, file_path, tmp_dir, stderr=sys.stderr, stdout=sys.stdout):
        self.file_path = file_path
        self.tmp_file = tempfile.mkstemp(dir=tmp_dir)[1]
        self.stderr = stderr
        self.stdout = stdout
        self.handle = None
        self.nucl = True

    def __enter__(self):
        self.handle = gzip.open(self.tmp_file, "wt")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.handle.close()
        self.index_seq()

    def write(self, record):
        if len(record.seq):
            bpio.write(record, self.handle, "fasta")
        else:
            self.stderr.write(f'{record.id} is empty')

    def index_seq(self):
        cmd = f'zcat "{self.tmp_file}" | bgzip > "{self.file_path}"'
        sp.run(cmd, shell=True, stderr=self.stderr, stdout=self.stdout)
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"{self.file_path} could not be created")
        os.remove(self.tmp_file)
        cmd = f'samtools faidx "{self.file_path}"'
        sp.run(cmd, shell=True, stderr=self.stderr, stdout=self.stdout)
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"{self.file_path[:-3] + '.fai'} could not be created")


class SeqQS:
    def __init__(self, db_file,stderr=sys.stderr):
        self.db_file = db_file
        self.seqs = SeqIO.index(db_file,"fasta")
        self.stderr = stderr

    def __getitem__(self, contig):
        return self.seq(contig)

    def __call__(self, contig, start=None, end=None):
        return self.seq(contig, start, end)

    def seq(self, contig, start=None, end=None, _num=0):
        region = contig
        if start and end:
            region = f'{contig}:{start}-{end}'
        elif start:
            region = f'{contig}:{start}'
        elif end:
            region = f'{contig}:-{end}'
        cmd = f'samtools faidx  "{self.db_file}" {region}'
        return "".join(sp.check_output(cmd, shell=True, stderr=self.stderr).decode("utf8").strip().split("\n")[1:])




class SeqStore:

    @staticmethod
    def instance():
        return SeqStore(SEQSTORE_DIR + ("" if SEQSTORE_DIR.endswith("/") else "/"),
                        SEQSTORE_DIR_TMP + ("" if SEQSTORE_DIR_TMP.endswith("/") else "/"),
                        )

    def __init__(self, base_path, temp_path, stderr=sys.stderr, stdout=sys.stdout):
        self.base_path = base_path
        self.temp_path = temp_path
        self.stderr = stderr
        self.stdout = stdout

    def db_path(self, accession):
        assert "/" not in accession, f"accession must not contain '/' : '{accession}' "
        return self.base_path + accession

    def genome_db_path(self, accession):
        return self.db_path(accession) + "/genome.fna.bgz"

    def proteome_db_path(self, accession):
        return self.db_path(accession) + "/proteome.faa.bgz"

    def genes_db_path(self, accession):
        return self.db_path(accession) + "/genes.fna.bgz"

    def add_genome(self, accession, record_iterator, force=False):
        self.add_seqs(accession, record_iterator, self.genome_db_path(accession), force)

    def add_seqs(self, accession, record_iterator, db: str, force=False,
                 stderr=sys.stderr, stdout=sys.stdout):

        if not force and os.path.exists(self.db_path(accession)):
            raise FileExistsError(
                f"{self.db_path(accession)} already exists, force the operation or delete the directory ")

        tmp_file = self.temp_path + accession + ".fasta"

        with self.stream(db, tmp_file, stdout=stdout, stderr=stderr) as stream:
            for record in record_iterator:
                stream.write(record)

    def stream(self, db, force=False,
               stderr=sys.stderr, stdout=sys.stdout):

        if not force and os.path.exists(db):
            raise FileExistsError(
                f"{db} already exists, force the operation or delete the directory ")

        assert db.endswith(".bgz"), f"{db} must end with .bgz"

        return SeqStream(db, self.temp_path, stderr=stderr, stdout=stdout)

    def make_blast(self, db, blast_name,
                   stderr=sys.stderr, stdout=sys.stderr):
        if not os.path.exists(db):
            raise FileNotFoundError(f"{db} does not exists")

        cmd = f'bgzip -d  -c {db} | makeblastdb -dbtype nucl -title {blast_name} -out {db} -in -'
        sp.run(cmd, shell=True, stderr=stderr, stdout=stdout)

    def merge_to(self, accession, general_db):
        """https://www.ncbi.nlm.nih.gov/books/NBK279693/
            blastdb_aliastool
        """
        raise NotImplemented()

    def qs(self, accession, seqtype: str = "genome"):
        assert seqtype in ["genome", "proteome"]
        db_file = self.genome_db_path(accession) if seqtype == "genome" else self.proteome_db_path(accession)
        return SeqQS(db_file)
