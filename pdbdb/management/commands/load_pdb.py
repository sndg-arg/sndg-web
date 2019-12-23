import os
import traceback
import warnings

from django.core.management.base import BaseCommand, CommandError

from pdbdb.models import PDB


from pdbdb.models import PDB
from pdbdb.io.PDBIO import PDBIO


class Command(BaseCommand):
    help = 'Imports a PDB'

    def add_arguments(self, parser):
        parser.add_argument('--code', required=True,help="4 letter PDB code")
        parser.add_argument('--tmp', default="/tmp/load_pdb")
        parser.add_argument('--pdbs_dir', default="/data/databases/pdb/divided/")
        parser.add_argument('--entries_path', default="/data/databases/pdb/entries.idx")

    def handle(self, *args, **options):
        pdbio = PDBIO(options['pdbs_dir'] + "/",options['entries_path'],options['tmp'])
        pdbio.init()
        try:
            pdbio.process_pdb(options['code'])

        except IOError as ex:
            traceback.print_exc()
            self.stderr.write("error processing pockets from %s: %s" % (pdb.code, str(ex)))
        except Exception as ex:
            traceback.print_exc()
            raise CommandError(ex)
