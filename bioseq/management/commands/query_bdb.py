import json
from json import JSONDecodeError

import Bio.SeqIO as bpio
from django.core.management.base import  CommandError
from tqdm import tqdm

from bioseq.models.Bioentry import Bioentry
from bioseq.models.Biodatabase import Biodatabase
from bioseq.io.SeqStore import SeqStore
from django_tqdm import BaseCommand

class Command(BaseCommand):
    help = 'Loads a genome in the database'

    def add_arguments(self, parser):
        parser.add_argument('--dbname', '-d', required=False)
        parser.add_argument('--query', '-q', default='{}')
        parser.add_argument('--format', '-f', default="auto", choices=["auto",'list', 'fasta'], )  # , 'gb', 'gff'
        parser.add_argument('--start', '-s', default=1, type=int)
        parser.add_argument('--end', '-e', default=None, type=int)

    def handle(self, *args, **options):
        dbname = options['dbname']
        query = options['query']
        out_format = options['format']

        if options['start'] < 1:
            raise CommandError(f'start must be greater than 0, "{options["start"]}" passed')
        if options['end'] and options['end'] < 0:
            raise CommandError(f'end must be greater than 0, "{options["start"]}" passed')

        options['start'] = options['start'] - 1
        
        
        
        if (out_format == "auto") and  dbname:
            out_format = "fasta"
        else:
            out_format = "list"

        if query:
            try:
                query = json.loads(query)
            except JSONDecodeError:
                if dbname:
                    query = {"name": query}
                else:
                    query = {"accession": query}

        if dbname:
            query["biodatabase__name"] = dbname
            query_manager = Bioentry.objects
        else:
            query_manager = Biodatabase.objects
        
        if not options['end']:
            if (out_format == "list") and (options['end'] == None):
                options['end'] = 10

        self.stderr.write(f"quering... {json.dumps(query)}")

        qs = query_manager.filter(**query)

        self.stderr.write(f"retreived sequences: {qs.count()}")

        seqstore = SeqStore.instance()
        seqtype = "genome"
        if dbname.endswith(Biodatabase.PROT_POSTFIX):
            dbname = dbname[:-len(Biodatabase.PROT_POSTFIX)]
            seqtype = "proteome"
        seq_qs = seqstore.qs(dbname,seqtype)
            
        qs2 = qs
        if options['start'] and options['end']:
            qs2 = qs2[options['start']:options['end']]
        elif options['start']:
            qs2 = qs2[options['start']:]
        elif options['end']:
            qs2 = qs2[:options['end']]

        for be in self.tqdm(qs2,  total=qs.count()):
            if out_format == "fasta":
                r = be.to_seq_record(seq_qs)
                bpio.write(r, self.stdout, out_format)
            else:
                self.stdout.write(be.name)

        if not dbname:
            if not options["end"]:
                options["end"] = qs.count()
            self.stderr.write(f'exported from {options["start"]} to {options["end"]} of {qs.count()}')

        self.stderr.write("finished!")
