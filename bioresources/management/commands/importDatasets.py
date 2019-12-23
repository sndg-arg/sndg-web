import json
import os
from datetime import datetime
from urllib.error import URLError

from Bio import Entrez
from django.core.management.base import BaseCommand
from django.db import transaction
from tqdm import tqdm

from bioresources.io.adapters import scopus_extended_publication, NCBIGDSAdapter, NCBIAssemblyAdapter, \
    NCBIBioSampleAdapter, NCBISRAAdapter, NCBIStructureAdapter
from bioresources.io.scopus import ScopusDS
from bioresources.models import ProcessStatus, ProcessStatusStep, Organization, ResourceRelation, Resource


def retry(q, n=3):
    for _ in range(n):
        try:
            data = q()
        except URLError as ex:
            continue
        return data
    raise ex


class Command(BaseCommand):
    """
    pip install -e git+https://github.com/ezequieljsosa/elsapy.git
    """

    help = 'Import resources from ncbi samples and publications for a given country'

    links = ["biosample", "assembly",
             "genome", "gds", "sra",
             "structure"]  # "bioproject" "nuccore"

    def add_arguments(self, parser):
        parser.add_argument('--scopus_config', required=True)
        parser.add_argument('--email', required=True)
        parser.add_argument('--country', required=True)
        parser.add_argument('--dst_dir', default="./")
        parser.add_argument('--job', default=datetime.now().strftime('%y-%m-%d'))
        parser.add_argument('--scopus_tmp', default="/tmp/scopus.json")

    def download_and_save_scopus(self, options):
        scopusDS = ScopusDS(options["scopus_config"])
        self.stdout.write("Downloading publications...")
        results = scopusDS.query(options["country"])
        first = True
        with open(options["scopus_tmp"], "w") as h:
            h.write("[")
            for x in results:
                if first:
                    first = False
                else:
                    h.write(",\n")
                json.dump(x, h)
            h.write("]")

    def cross_pubmed_ncbi_data(self, ncbi_id):

        pmc = Entrez.read(retry(lambda: Entrez.esummary(db="pubmed", id=ncbi_id)))[0]["ArticleIds"]

        if "pmc" in pmc:
            pmc_id = pmc["pmc"].replace("PMC", "")
            # https://www.ncbi.nlm.nih.gov/pmc/?Db=nuccore&DbFrom=pmc&Cmd=Link&LinkName=pmc_nuccore&IdsFromResult=5814494
            with retry(lambda: Entrez.elink(dbfrom="pmc", id=pmc_id, db=",".join(Command.links))) as handle:
                links = Entrez.read(handle)[0]
        else:
            with retry(lambda: Entrez.elink(dbfrom="pubmed", id=str(ncbi_id), db=",".join(Command.links))) as handle:
                links = Entrez.read(handle)[0]

        data = []
        if links["LinkSetDb"]:
            for linklist in links["LinkSetDb"]:
                ids = [str(list(x.values())[0]) for x in linklist["Link"]]
                record = {"ncbi_publication_id": ncbi_id,
                          "ncbi_db_type": str(linklist["LinkName"]).split("_")[1],
                          "ncbi_db_id": ids}
                data.append(record)
        return data

    def cross_biosample_ncbi_data(self, ncbi_id):

        with retry(lambda: Entrez.elink(dbfrom="biosample", id=ncbi_id, db=",".join(Command.links))) as handle:
            pmc_record = Entrez.read(handle)[0]

        data = []
        for linklist in pmc_record["LinkSetDb"]:
            ids = [str(list(x.values())[0]) for x in linklist["Link"]]
            record = {"ncbi_sample_id": ncbi_id, "ncbi_db_type": linklist["LinkName"].split("_")[1],
                      "ncbi_db_id": ids}
            data.append(record)
        return data




    def handle(self, *args, **options):

        self.mapping = {
            "gds": NCBIGDSAdapter,
            "structure": NCBIStructureAdapter,
            "sra": NCBISRAAdapter,
            "assembly": NCBIAssemblyAdapter,
            "biosample": NCBIBioSampleAdapter,

        }

        assert os.path.exists(options["dst_dir"])
        assert os.path.exists(options["scopus_config"])

        options["scopus_config"] = os.path.abspath(options["scopus_config"])
        options["dst_dir"] = os.path.abspath(options["dst_dir"])
        Entrez.email = options["email"]
        Organization.objects.get_or_create(name="NCBI")

        ps = ProcessStatus.objects.filter(name=options["job"]).first()
        if not ps:
            ps = ProcessStatus(name=options["job"])
            ps.save()
            step = ProcessStatusStep(name="load_publications", process_status=ps,
                                     class_identifier="id", class_name="bioresources.models.Publication")
            step.save()

            step = ProcessStatusStep(name="processed_publication", process_status=ps,
                                     class_identifier="id", class_name="bioresources.models.Publication")
            step.save()

            step = ProcessStatusStep(name="country_samples", process_status=ps,
                                     class_identifier="id", class_name="bioresources.models.Sample")
            step.save()

        if not os.path.exists(options["scopus_tmp"]):
            self.download_and_save_scopus(options)
        load_publications_step = ps.step("load_publications")
        if not ps.step("load_publications").completed:
            with tqdm(json.load(open(options["scopus_tmp"]))) as pbar:
                added_publications = 0

                for article in pbar:
                    pbar.set_description("Processing Publications. Added %i of " % added_publications)
                    publication_id = None
                    if article['dc:identifier'] not in load_publications_step:

                        with transaction.atomic():
                            publication = scopus_extended_publication(article)
                            if publication:
                                publication_id = publication.id
                                added_publications += 1
                    load_publications_step.append(publication_id, article['dc:identifier'])

            load_publications_step.completed = True
            load_publications_step.save()

        processed_publication_step = ps.step("processed_publication")
        if not processed_publication_step.completed:
            with tqdm(load_publications_step.results(), total=load_publications_step.units.count()) as pbar:
                pbar.set_description("Fetching paper datasets from NCBI")
                for publication in pbar:
                    if publication.id not in processed_publication_step:
                        if publication.doi:
                            with retry(lambda: Entrez.esearch(db="pubmed", retmax=10, term=publication.doi)) as h:
                                pubmed_record = Entrez.read(h)
                                if pubmed_record["IdList"]:
                                    ncbi_id = pubmed_record["IdList"][0]
                                    crossed_data_sets = self.cross_pubmed_ncbi_data(ncbi_id)
                                    for crossed_data in crossed_data_sets:
                                        with transaction.atomic():
                                            for crossed_id in crossed_data["ncbi_db_id"]:
                                                if crossed_data["ncbi_db_type"] in self.mapping:
                                                    mapper = self.mapping[crossed_data["ncbi_db_type"]]()
                                                    data = mapper.fetch(crossed_id)
                                                    dataset = mapper.save(data, crossed_id)
                                                    (ResourceRelation.
                                                        objects.get_or_create(source=publication,
                                                                              target=dataset,
                                                                              role=Resource.RESOURCE_TYPES[
                                                                                       Resource.RESOURCE_TYPES.PUBLICATION] +
                                                                                   "_" + Resource.RESOURCE_TYPES[
                                                                                       dataset.type]))
                        processed_publication_step.append(publication.id)
            processed_publication_step.completed = True
            processed_publication_step.save()

        # By country sample
        country_samples_step = ps.step("country_samples")
        if not country_samples_step.completed:
            search = '{country}[All Fields] AND "attribute geographic location"[filter]'.format(
                country=options["country"])
            esearch = Entrez.read(retry(lambda: Entrez.esearch(db="biosample", term=search, retmax=10000)))["IdList"]
            for biosample_id in tqdm(esearch):
                biosample_id = str(biosample_id)
                if biosample_id not in country_samples_step:
                    dataset_id = None
                    with transaction.atomic():
                        links = Entrez.read(
                            retry(
                                lambda: Entrez.elink(dbfrom="biosample", id=biosample_id, db=",".join(Command.links))))

                        if len(links[0]['LinkSetDb']) > 0:
                            mapper = self.mapping["biosample"]()
                            biosampledata = mapper.fetch(biosample_id)
                            records = biosampledata["Attributes"]["Attribute"]
                            if isinstance(records, dict):
                                records = [records]
                            attributes = {x["@attribute_name"]: x["#text"] for x in records}

                            sample_geo_props = ['geographic location (country and/or sea)', 'geo_loc_name',
                                                'Origin (developed or donated from)', "birth_location",
                                                "geo-loc-name", 'country_of_birth']

                            if any([options["country"].lower() in attributes.get(x, "").lower() for x in
                                    sample_geo_props]):
                                biosample = mapper.save(biosampledata, biosample_id)
                                dataset_id = biosample.id
                                crossed_data_sets = self.cross_biosample_ncbi_data(biosample_id)
                                for crossed_data in crossed_data_sets:
                                    for crossed_id in crossed_data["ncbi_db_id"]:
                                        mapper = self.mapping[crossed_data["ncbi_db_type"]]()
                                        data = mapper.fetch(crossed_id)
                                        dataset = mapper.save(data, crossed_id)
                                        ResourceRelation.objects.get_or_create(source=biosample, target=dataset,
                                                                               role=
                                                                               Resource.RESOURCE_TYPES[
                                                                                   Resource.RESOURCE_TYPES.SAMPLE] +
                                                                               "_" + Resource.RESOURCE_TYPES[
                                                                                   dataset.type])
                            else:
                                print([k for k, v in attributes.items() if options["country"] in v.lower()])
                    country_samples_step.append(str(dataset_id), biosample_id)

            country_samples_step.completed = True
            country_samples_step.save()
