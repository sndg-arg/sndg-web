import os
import sys
import re
import subprocess as sp
import warnings

from Bio import BiopythonWarning, BiopythonParserWarning, BiopythonDeprecationWarning, BiopythonExperimentalWarning
from django.core.management.base import BaseCommand, CommandError

from bioseq.io.BioIO import BioIO
from bioseq.io.SeqStore import SeqStore

warnings.simplefilter('ignore', RuntimeWarning)
warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)
warnings.simplefilter('ignore', BiopythonDeprecationWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

from SNDG.Sequence import smart_parse
from SNDG import mkdir



class Command(BaseCommand):
    help = 'Loads a genome in the database'

    def add_arguments(self, parser):
        parser.add_argument('--input', '-i', required=True)
        parser.add_argument('--accession', '-a', required=True)
        parser.add_argument('--description', '-d', required=False)
        parser.add_argument('--seqs', '-s', default=None)
        parser.add_argument('--taxon', '-t', type=int, required=False)
        parser.add_argument('--force', action="store_true")

    def handle(self, *args, **options):
        input_file = options['input']
        accession = options['accession']

        self.stderr.write(f"trying to import  {options['accession']} imported!")

        if not os.path.exists(input_file):
            raise CommandError(f'{input_file} does not exists')

        extra_attrs = {}
        taxon = self.detect_tax(input_file, extra_attrs)
        taxon = options['taxon'] if options['taxon'] else taxon
        description = options["description"] if options["description"] else " ".join(
            [f'{k}:{v}' for k, v in extra_attrs.items()])

        io = BioIO(accession, taxon, stderr=self.stderr)

        if options["force"]:
            if io.exists():
                res = io.delete()
                self.stderr.write(str(res))
            elif io.exists():
                raise CommandError(f'{accession} already exists, use --force to overwrite ')

        grep_cmd = 'grep -c "FEATURES *Location/Qualifiers" "%s"' % input_file
        if input_file.endswith(".gz"):
            grep_cmd = 'z' + grep_cmd
        total = int(sp.check_output(grep_cmd, shell=True))
        io.create_db(description)

        seqstore = SeqStore.instance()

        if options['seqs']:
            if not os.path.exists(options['seqs']):
                raise CommandError(f'{options["seqs"]} does not exists')

            it = smart_parse(input_file, smart_parse(options["seqs"]))
        else:
            it = smart_parse(input_file)

        mkdir(seqstore.db_path(accession))
        s1 = seqstore.stream(seqstore.genome_db_path(accession), force=True,
                             stderr=sys.stderr, stdout=sys.stderr)
        s2 = seqstore.stream(seqstore.proteome_db_path(accession), force=True,
                             stderr=sys.stderr, stdout=sys.stderr)
        with s1 as genome_stream, s2 as proteome_stream:
            io.process_record_list(it, total, genome_stream, proteome_stream)

        self.stderr.write(f"genome {options['accession']} imported!")

    def detect_tax(self, input_file, other_attrs: dict):
        # /strain="AYE"
        # /plasmid="
        # /organism="
        # ;old-name=Acinetobacter baumannii;plasmid-name=p1ABAYE;strain=AYE
        if re.sub(r'.gz$', '', input_file).split(".")[-1] in ["gb", "gbk", "gbff"]:
            # /db_xref="taxon:509173"
            search_txt = '/db_xref="taxon:'
            search_txt2 = '/strain="'
            search_txt3 = '/organism="'
        elif input_file.split(".")[-1].replace(".gz", "") in ["gff", "gtf"]:
            # ID=NC_010410.1:1..3936291;Dbxref=taxon:509173;Is_circul
            search_txt = "Dbxref=taxon:"
            search_txt2 = 'strain='
            search_txt3 = 'old-name='
        else:
            raise CommandError(f'could not detect tax in {input_file}')
        cat = ("z" if input_file.endswith(".gz") else "") + "cat"
        grep_cmd = f'''{cat} "{input_file}" | head -n 100 | grep  '{search_txt}' '''
        try:
            tax_line = sp.check_output(grep_cmd, shell=True).decode("utf-8").strip().split("\n")[0]
            tax = tax_line.split(search_txt)[1].split(";")[0].replace('"', "")
            tax = int(tax)
        except:
            raise CommandError(f'could not detect tax in {input_file}')

        try:
            grep_cmd = f'''{cat} "{input_file}" | head -n 100 | grep  '{search_txt2}' '''
            line = sp.check_output(grep_cmd, shell=True).decode("utf-8").strip().split("\n")[0]
            other_attrs["strain"] = line.split(search_txt2)[1].split(";")[0].replace('"', "")

        except:
            pass

        try:
            grep_cmd = f'''{cat} "{input_file}" | head -n 100 | grep  '{search_txt3}' '''
            line = sp.check_output(grep_cmd, shell=True).decode("utf-8").strip().split("\n")[0]
            other_attrs["organism"] = line.split(search_txt3)[1].split(";")[0].replace('"', "")

        except:
            pass

        return tax
