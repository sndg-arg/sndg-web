import os

from django.core.management.base import BaseCommand, CommandError
from tqdm import tqdm

from pdbdb.io.PDB2SQL import PDB2SQL
from pdbdb.models import PDB


def iterpdbs(pdbs_dir, pdb_extention=".ent"):
    for index_dir in os.listdir(pdbs_dir):
        if len(index_dir) == 2:
            for x in os.listdir(pdbs_dir + "/" + index_dir):
                if x.endswith(pdb_extention):
                    yield x[3:7], pdbs_dir + "/" + index_dir + "/" + x


class Command(BaseCommand):
    help = 'Loads the pdb files to the database'

    def add_arguments(self, parser):
        parser.add_argument('--pdbs_dir', required=True)
        parser.add_argument('--entries_path', required=True)

    def handle(self, *args, **options):

        pdb2sql = PDB2SQL(options['pdbs_dir'], options['entries_path'])
        pdb2sql.load_entries()
        pdbs = list(tqdm(iterpdbs(pdbs_dir)))
        # 4zux 42 mer 2lo7("5my5","/data/databases/pdb/divided/my/pdb5my5.ent")
        # ("4zu4", "/data/databases/pdb/divided/zu/pdb4zu4.ent")

        with tqdm(pdbs) as pbar:
            for code in pbar:
                code = code.lower()
                pdb_path = pdb2sql.download(code)

                if PDB.objects.filter(code=code).exists():
                    continue

                pbar.set_description(code)
                try:
                    pdb2sql.create_pdb_entry(code, pdb_path)
                    pdb2sql.update_entry_data(code, pdb_path)
                except Exception as ex:
                    raise CommandError(ex)
