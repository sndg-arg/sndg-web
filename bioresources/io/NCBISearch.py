import os
import subprocess as sp

from Bio import Entrez

from .adapters import NCBIAssemblyAdapter, NCBIPubmedAdapter, NCBIPmcAdapter, NCBIBioSampleAdapter, NCBIBioProject, \
    NCBIStructureAdapter, NCBIGDSAdapter, NCBISRAAdapter, NCBINuccoreAdapter, NCBIGeneAdapter, NCBIProteinAdapter

from bioresources.models.Resource import Resource

from django.core.cache import cache


class NCBISearch():
    allowed_dbs = ["gds", "pmc", "pubmed", "bioproject", "biosample", "sra", "taxonomy", "structure", "genome",
                   "protein", "gene", "nuccore", "assembly"]
    database_map = {"assembly": NCBIAssemblyAdapter(), "pmc": NCBIPmcAdapter(), "pubmed": NCBIPubmedAdapter(),
                    "gene": NCBIGeneAdapter(), "protein": NCBIProteinAdapter(), "structure": NCBIStructureAdapter(),
                    "gds": NCBIGDSAdapter(), "bioproject": NCBIBioProject(),
                    "nuccore": NCBINuccoreAdapter(), "sra": NCBISRAAdapter(), "biosample": NCBIBioSampleAdapter()

                    }
    db_type = {
        "assembly": Resource.RESOURCE_TYPES.ASSEMBLY, "pmc": Resource.RESOURCE_TYPES.PUBLICATION,
        "pubmed": Resource.RESOURCE_TYPES.PUBLICATION,
        "gene": 40, "protein": 40, "structure": Resource.RESOURCE_TYPES.STRUCTURE,
        "gds": Resource.RESOURCE_TYPES.EXPRESSION, "bioproject": Resource.RESOURCE_TYPES.BIOPROJECT,
        "nuccore": 40, "sra": Resource.RESOURCE_TYPES.READS, "biosample": Resource.RESOURCE_TYPES.SAMPLE
    }

    rtype2ncbb = {v: k for k, v in db_type.items()}

    def search_all(self, search):
        result = {}
        data = Entrez.read(Entrez.egquery(term=search))["eGQueryResult"]
        for dbResult in data:
            if dbResult["DbName"] in NCBISearch.allowed_dbs:
                try:
                    result[dbResult["DbName"]] = int(dbResult["Count"])
                except:
                    result["DbName"] = 0
        result["assembly"] = len(
            Entrez.read(Entrez.esearch(db="assembly", term=search + "[All Names] OR " + search + "[All Uids]"))[
                'IdList'])
        return result


    def search_database(self, database, search, retmax=20, retstart=0):
        records = {}
        adapter = NCBISearch.database_map[database]
        esearch = Entrez.read(Entrez.esearch(db=database, term=search, retmax=retmax, retstart=retstart))
        id_list = list(esearch["IdList"])
        if id_list:
            for ncbi_id, data in zip(id_list, adapter.fetch_list(",".join(id_list))):
                records[str(ncbi_id)] = adapter.adapt(data)
            return records, int(esearch["Count"])
        else:
            return {}, 0

    def save_resource(self, database: str, ncbi_id: str) -> Resource:

        adapter = NCBISearch.database_map[database]
        data = adapter.fetch(ncbi_id)
        return adapter.save(data, ncbi_id)


    @staticmethod
    def download_assembly( assembly_name, workdir):
        arr = assembly_name.split("_")
        acc = "_".join(arr[:2])
        name = "_".join(arr[2:])
        asspath = "/".join([acc[0:3], acc[4:7], acc[7:10], acc[10:13],
                            acc + "_" + name.replace(" ", "_").replace("#", "_")])

        cmd = 'rsync --recursive --include="*_genomic.gbff.gz" --exclude="*"   rsync://ftp.ncbi.nlm.nih.gov/genomes/all/' + asspath + '/ "' + \
              workdir + '"'
        with open(os.devnull, 'w') as FNULL:
            sp.call(cmd, shell=True,stdout=FNULL)

    @staticmethod
    def doi(doi):
        result = Entrez.read(Entrez.esearch("pmc",doi + '[DOI]'))
        if "IdList" in result:
            ncbi_id = result["IdList"][0]
            publ = Entrez.read(Entrez.esummary(db="pmc",id=ncbi_id))
        else:
            return None