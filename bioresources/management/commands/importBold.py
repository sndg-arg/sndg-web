import json
import os
import subprocess as sp
import sys

from bioseq.models.Taxon import Taxon
from django.core.management.base import BaseCommand
from django.db import transaction
from tqdm import tqdm

from bioresources.models.Barcode import Barcode
from bioresources.models.Organization import Organization
from bioseq.models.Biodatabase import Biodatabase
from bioseq.models.Bioentry import Bioentry
from bioseq.models.Biosequence import Biosequence


def execute(cmd, **kwargs):
    sp.call(cmd.format(**kwargs), shell=True)


def download_file(complete_url, target, ovewrite=False, retries=3):
    if not target.strip():
        target = "./"
    if not os.path.exists(os.path.dirname(os.path.abspath(target))):
        raise Exception("%s does not exists" % os.path.dirname(target))
    if os.path.exists(target) and not ovewrite:
        raise OvewriteFileException("%s already exists" % target)

    execute(' wget  --timeout=20 --tries={retries} -O {target} "{url}"',
            url=complete_url, retries=retries, target=target)


class Command(BaseCommand):
    help = 'Download and loads all barcodes of Bold of a given country'

    def __init__(self, stdout=None, stderr=None, no_color=False):
        super().__init__(stdout=stdout, stderr=stderr, no_color=no_color)
        self.ontology = None
        self.dbmap = {}
        self.is_a = None
        self.relmap = {}
        self.ontology_name = ""
        self.obo_path = None

    def add_arguments(self, parser):

        parser.add_argument('--country', required=True)
        parser.add_argument('--bold_url',
                            default="http://www.boldsystems.org/index.php/API_Public/combined?format=json&geo=")
        parser.add_argument('--json', default="/tmp/bold.json")
        parser.add_argument('--override', action='store_false', help="override json tmp file. Default: true")
        parser.add_argument('--resources', action='store_false', help="update Barcodes. Default: true")
        parser.add_argument('--seqs', action='store_false', help="update sequences. Default: true")

    def handle(self, *args, **options):

        if options ["override"] or (not os.path.exists(options["json"])):
            download_file(options["bold_url"] + options["country"], options["json"],ovewrite=options ["override"])
        if options["resources"]:
            self.update_resources(options["json"])
        if options["seqs"]:
            self.update_biodb(options["json"])

    def update_biodb(self, json_input):
        bdb = Biodatabase.objects.get_or_create(name=Barcode.BIODBNAME)[0]
        with open(json_input) as h:
            data = json.load(h)["bold_records"]["records"].values()
        for d in tqdm(data):
            if "sequences" in d and Barcode.objects.filter(name=d["processid"]).count():
                if not Bioentry.objects.filter(biodatabase=bdb,accession=d["processid"]).count():
                    with transaction.atomic():
                        be = Bioentry(biodatabase=bdb, accession=d["processid"], identifier=d["processid"],
                                      name=d["processid"])
                        be.save()
                        seq = d["sequences"]["sequence"][0]["nucleotides"].replace("-", "")
                        Biosequence(bioentry=be, seq=seq, length=len(seq)).save()

    def update_resources(self, json_input):
        with open(json_input) as h:
            data = json.load(h)["bold_records"]["records"].values()

        bcodes = []
        org_bold = Organization.objects.get_or_create(name="BOLD")[0]

        for d in tqdm(data):
            if ("sequences" in d) and (not Barcode.objects.filter(name=d["processid"]).count()):

                for x in ["subspecies","species", "genus", "subfamily", "family", "order", "class", "phylum"]:
                    if x in d["taxonomy"]:
                        tax = d["taxonomy"][x]["taxon"]["taxID"]
                        break

                d["description"] = (d["taxonomy"][x]["taxon"]["name"] + " " +
                                    d["sequences"]["sequence"][0]["markercode"] + " "
                                    + d["specimen_identifiers"]["institution_storing"])

                d["tax"] = int(tax)

                bc = Barcode(
                    name=d["processid"],
                    description=d["description"],
                    country=d["collection_event"]["country"],
                    collectors=d["collection_event"]["collectors"],
                    type="barcode",
                    marker=d["sequences"]["sequence"][0]["markercode"],

                    bold_org=d["specimen_identifiers"]["institution_storing"]
                )

                for x in ["subspecies","species", "genus", "subfamily", "family", "order", "class", "phylum"]:
                    if x in d["taxonomy"]:
                        try:
                            tax = d["taxonomy"][x]["taxon"]["taxID"]
                            bc.ncbi_tax = Taxon.objects.get(ncbi_taxon_id=int(tax))
                        except Taxon.DoesNotExist:
                            continue
                        break


                try:
                    bc.image_url = d["specimen_imagery"]["media"][0]["image_file"]
                except KeyError:
                    pass
                try:
                    bc.subdivision = d["collection_event"]["province_state"]
                except KeyError:
                    pass

                bcodes.append(bc)
        bcodes_col = []
        for i, bc in enumerate(tqdm(bcodes)):
            bcodes_col.append(bc)
            if i == 5000:
                with transaction.atomic():
                    for bc2 in bcodes_col:
                        bc2.save()
                        bc2.publishers.add(org_bold)
                    bcodes_col = []

        with transaction.atomic():
            for bc2 in bcodes_col:
                bc2.save()
                bc2.publishers.add(org_bold)
