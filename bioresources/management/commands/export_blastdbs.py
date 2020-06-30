import os
from django.conf import settings

from tqdm import tqdm
import Bio.SeqIO as bpio

from django.core.management.base import BaseCommand, CommandError

from bioseq.io.SeqStore import SeqStore
from bioseq.models.Bioentry import Bioentry
from bioresources.models.BlastDB import BlastDB
from bioseq.io.Blast import Blast


class Command(BaseCommand):
    help = 'exports fasta files'

    def add_arguments(self, parser):
        parser.add_argument('--blastdbs_dir', '-p', default=getattr(settings, "BLASTDBSDIR"))
        for db in BlastDB.dbs:
            parser.add_argument('--' + db, action='store_false')

    def handle(self, *args, **options):
        seqstore = SeqStore.instance()
        for db, db_filter in BlastDB.dbs.items():

            if options[db]:
                bdb = BlastDB.objects.get(accession=db)
                file_path = options["blastdbs_dir"] + "/%s.fasta" % db
                with open(file_path, "w") as h:
                    qs = Bioentry.objects.filter(db_filter)
                    for be in tqdm(qs.all(), total=qs.count(), file=self.stderr):
                        bpio.write(be.to_seq_record(seqstore,False), h, "fasta")
                Blast.makeblastdb(fasta_path=file_path, dbtype=bdb.dbtype)
                bdb.save()
