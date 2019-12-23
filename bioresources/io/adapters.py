import datetime
from urllib.error import URLError
import xmltodict

from django.db import transaction

from Bio import Entrez
from bioseq.models.Taxon import Taxon, TaxonName
from bioseq.models.Term import Term
from bioseq.models.Ontology import Ontology
from bioseq.models.Bioentry import Bioentry

from bioresources.models.Assembly import Assembly
from bioresources.models.Sample import Sample
from bioresources.models.Publication import Publication
from bioresources.models.ExternalId import ExternalId
from bioresources.models.Resource import Resource
from bioresources.models.Organization import Organization
from bioresources.models.ReadsArchive import ReadsArchive
from bioresources.models.Structure import Structure
from bioresources.models.Expression import Expression
from bioresources.models.Affiliation import Affiliation
from bioresources.models.BioProject import BioProject
from bioresources.models.Person import Person

from bioresources.models.ResourceProperty import ResourceProperty, ResourcePropertyValue


def retry(q, n=4):
    for _ in range(n):
        try:
            data = q()
        except URLError as ex:
            continue
        return data
    raise ex


def ncbi_link(dbfrom, id, linkname):
    pass


def scopus_publication(article):
    if Publication.objects.filter(scopus_id=article['dc:identifier']).exists():
        publication = Publication.objects.get(scopus_id=article['dc:identifier'])
    elif ("prism:doi" in article) and \
            Publication.objects.filter(doi=article['prism:doi']).exists():
        publication = Publication.objects.get(doi=article['prism:doi'])
    elif Publication.objects.filter(name=article["dc:title"][:350]).exists():
        publication = Publication.objects.get(name=article["dc:title"][:350])
    else:
        publication = Publication(
            type=Resource.RESOURCE_TYPES.PUBLICATION,
            name=article["dc:title"][:350],
            date_of_publication=datetime.datetime.strptime(article['prism:coverDate'], "%Y-%m-%d"),
        )
        publication.scopus_id = article['dc:identifier']
        publication.electronic_id = article["eid"]
        if "dc:description" in article:
            publication.description = article["dc:description"]
        if 'pubmed-id' in article:
            publication.pubmed_id = article['pubmed-id']
        if "prism:doi" in article:
            publication.doi = article["prism:doi"]
        if 'prism:issn' in article:
            publication.issn = article['prism:issn']

    return publication


def scopus_affiliation(affiliation):
    afcountry = affiliation["affiliation-country"]

    org = Organization(name=affiliation["affilname"],
                       country=afcountry,
                       city=affiliation["affiliation-city"],
                       source=Organization.objects.get(name=Organization.SCOPUS)
                       )

    # TODO country detection based on the city
    if not affiliation["affiliation-country"]:
        if any([x in str(org.city) for x in [
            "CABA", "Buenos Aires", "Rosario"
        ]]):
            org.country = "Argentina"

    if "afid" in affiliation:
        org.scopus_id = affiliation['afid']
        if Organization.objects.filter(scopus_id=org.scopus_id).exists():
            org = Organization.objects.get(scopus_id=org.scopus_id)
        else:
            org.save()

    return org


def scopus_author(author, publication, arg):
    person = Person(surname=author["surname"],
                    name=author["given-name"] if author["given-name"] else "",
                    scopus_id=author["authid"])
    if Person.objects.filter(scopus_id=person.scopus_id).exists():
        person = Person.objects.get(scopus_id=person.scopus_id)
    else:
        person.save()

    if ("afid" in author):

        aff = Affiliation(resource=publication, author=person)
        aff.save()
        for affdict in author["afid"]:
            aff.organizations.add(Organization.objects.get(scopus_id=affdict["$"]))
        aff.save()

        if [x for x in author["afid"] if x["$"] in arg]:
            person.save()
    return person


def scopus_extended_publication(article):
    publication = scopus_publication(article)
    publication.save(force_insert=True)
    arg = []  # TODO criterio pais?
    for affiliation in article["affiliation"]:
        org = scopus_affiliation(affiliation)
        #     if org.country == "Argentina" and org.scopus_id:
        arg.append(org.scopus_id)
    #
    # if not arg:
    #     # for a in article["affiliation"]:
    #     #     if not a["affiliation-country"]:
    #     #         if a["affiliation-city"]:
    #     #             pepe.append(a["affiliation-city"])
    #     return

    if "author" in article:
        for author in article["author"]:
            scopus_author(author, publication, arg)
    return publication


class NCBIBioSampleAdapter:

    def fetch(self, ncbi_id):
        return self.fetch_list(ncbi_id)[0]

    def fetch_list(self, ncbi_ids):
        slist = Entrez.read(retry(lambda: Entrez.esummary(db="biosample", id=ncbi_ids)),validate=False)["DocumentSummarySet"][
            "DocumentSummary"]
        return [xmltodict.parse(biosampledata["SampleData"])["BioSample"] for biosampledata in slist]

    def adapt(self, summaryData) -> Sample:

        acc = summaryData["@accession"]
        desc = summaryData["Description"]["Title"]

        s = Sample(name=acc, description=desc, type=Sample.TYPE)

        if ("Organism" in summaryData["Description"]) and ("@taxonomy_id" in summaryData["Description"]["Organism"]):
            tax = summaryData["Description"]["Organism"]['@taxonomy_id']
            s.ncbi_tax = Taxon.objects.filter(ncbi_taxon_id=int(tax)).first()

        try:
            s.publication_date = datetime.datetime.strptime(summaryData["@publication_date"].split("T")[0], "%Y-%m-%d")
        except ValueError:
            pass
        try:
            s.update_date = datetime.datetime.strptime(summaryData["@last_update"].split("T")[0], "%Y-%m-%d")
        except ValueError:
            pass
        return s

    def save(self, summaryData, ncbi_id):

        s = self.adapt(summaryData)
        s.save()
        ncbi_org = Organization.objects.get(name="NCBI")
        s.publishers.add(ncbi_org)

        self.save_attributes(summaryData,s)

        ExternalId(resource=s, organization=ncbi_org,
                   identifier=s.name, type="accession").save()
        ExternalId(resource=s, organization=ncbi_org,
                   identifier=ncbi_id, type="identifier").save()
        return s

    def save_attributes(self, summaryData,s):
        ncbi_org = Organization.objects.get(name="NCBI")
        ontology = Ontology.objects.get_or_create(name="NCBI sample", definition="Attributes of an NCBI Sample")[0]
        records = summaryData["Attributes"]["Attribute"]
        if isinstance(records, dict):
            records = [records]

        for x in records:
            if x["#text"].strip() and (x["#text"] != "."):
                term = Term.objects.get_or_create(
                    ontology=ontology, name=x["@attribute_name"], identifier=x["@attribute_name"][:255])[0]
                prop = ResourceProperty.objects.get_or_create(term=term, resource=s, organization=ncbi_org)[0]
                ResourcePropertyValue.objects.create(property=prop, value=x["#text"][:200])


class NCBISRAAdapter:

    def fetch(self, ncbi_id):
        return Entrez.read(retry(lambda: Entrez.esummary(db="sra", id=ncbi_id)))[0]

    def fetch_list(self, ncbi_ids):
        return Entrez.read(retry(lambda: Entrez.esummary(db="sra", id=ncbi_ids)))

    def adapt(self, summaryData):

        # sra_record = Entrez.read(retry(lambda: Entrez.esummary(db="sra", id=ncbi_id)))
        expData = xmltodict.parse("<xml>" + summaryData["ExpXml"] + "</xml>")["xml"]
        acc = expData["Experiment"]["@acc"]

        qs = ReadsArchive.objects.filter(external_ids__identifier=acc, external_ids__type="accession")
        if qs.exists():
            return qs.first()

        s = ReadsArchive(type=Resource.RESOURCE_TYPES.READS, name=acc, description=expData["Summary"]["Title"])
        s.id = Resource.objects.latest('id').id + 1
        try:
            s.release_date = datetime.datetime.strptime(summaryData["CreateDate"].split(" ")[0], "%Y/%m/%d")
        except ValueError:
            pass
        try:
            s.update_date = datetime.datetime.strptime(summaryData["UpdateDate"].split(" ")[0], "%Y/%m/%d")
        except ValueError:
            pass

        if ("Organism" in expData) and ("@taxid" in expData["Organism"]):
            s.ncbi_tax = Taxon.objects.filter(ncbi_taxon_id=int(expData["Organism"]["@taxid"])).first()
        return s

    def save(self, summaryData, ncbi_id):
        s = self.adapt(summaryData)

        ncbi_org = Organization.objects.get(name="NCBI")
        s.type = s.__class__.TYPE
        s.save(force_insert=True)
        s.publishers.add(ncbi_org)

        ExternalId(resource=s, organization=ncbi_org,
                   identifier=s.name, type="accession").save(force_insert=True)
        ExternalId(resource=s, organization=ncbi_org,
                   identifier=ncbi_id, type="identifier").save(force_insert=True)
        return s


class NCBIPmcAdapter:
    def fetch(self, ncbi_id):
        return self.fetch_list(ncbi_id)[0]

    def fetch_list(self, ncbi_ids):
        return Entrez.read(retry(lambda: Entrez.esummary(db="pmc", id=ncbi_ids)))

    def adapt(self, summaryData, ncbi_id):
        s = Publication(name=summaryData["Title"], description="")
        s.type = s.__class__.TYPE
        return s


class NCBIPubmedAdapter:
    def fetch(self, ncbi_id):
        return self.fetch_list(ncbi_id)[0]

    def fetch_list(self, ncbi_ids):
        return Entrez.read(retry(lambda: Entrez.esummary(db="pubmed", id=ncbi_ids)))

    def adapt(self, summaryData, ncbi_id):
        s = Publication(name=summaryData["Title"], description="")
        s.type = s.__class__.TYPE
        return s


class NCBIGeneAdapter:
    def fetch(self, ncbi_id):
        return self.fetch_list(ncbi_id)[0]

    def fetch_list(self, ncbi_ids):
        return Entrez.read(retry(lambda: Entrez.esummary(db="gene", id=ncbi_ids)))["DocumentSummarySet"][
            "DocumentSummary"]

    def adapt(self, summaryData, ncbi_id):
        s = Bioentry(name=summaryData["Name"], description=summaryData["Summary"])
        return s


class NCBIProteinAdapter:
    def fetch(self, ncbi_id):
        return self.fetch_list(ncbi_id)[0]

    def fetch_list(self, ncbi_ids):
        return Entrez.read(retry(lambda: Entrez.esummary(db="protein", id=ncbi_ids)))["DocumentSummarySet"][
            "DocumentSummary"]

    def adapt(self, summaryData):
        return Bioentry(name=summaryData["title"], description=summaryData["abstract"])


class NCBINuccoreAdapter:
    def fetch(self, ncbi_id):
        return self.fetch_list(ncbi_id)[0]

    def fetch_list(self, ncbi_ids):
        return Entrez.read(retry(lambda: Entrez.esummary(db="nuccore", id=ncbi_ids)))

    def adapt(self, summaryData):
        return Bioentry(name=str(summaryData["Gi"]), description=summaryData["Title"])


class NCBIBioProject:
    def fetch(self, ncbi_id):
        return self.fetch_list(ncbi_id)[0]

    def fetch_list(self, ncbi_ids):
        return Entrez.read(retry(lambda: Entrez.esummary(db="bioproject", id=ncbi_ids)))["DocumentSummarySet"][
            "DocumentSummary"]

    def adapt(self, summaryData):
        s = BioProject(name=str(summaryData["Project_Acc"]), description=summaryData["Project_Title"])
        s.type = s.__class__.TYPE
        return s

    def save(self, summaryData, ncbi_id):
        s = self.adapt(summaryData)
        s.save(force_insert=True)

        ncbi_org = Organization.objects.get(name="NCBI")
        s.publishers.add(ncbi_org)

        ExternalId(resource=s, organization=ncbi_org,
                   identifier=s.name, type="accession").save(force_insert=True)
        ExternalId(resource=s, organization=ncbi_org,
                   identifier=ncbi_id, type="identifier").save(force_insert=True)
        return s


class NCBIAssemblyAdapter:

    def fetch(self, ncbi_id):
        return self.fetch_list(ncbi_id)[0]

    def fetch_list(self, ncbi_ids):
        return Entrez.read(retry(lambda: Entrez.esummary(db="assembly", id=ncbi_ids)))["DocumentSummarySet"][
            "DocumentSummary"]

    def adapt(self, summaryData):
        name = summaryData["AssemblyName"]
        acc = summaryData["AssemblyAccession"]

        tax = Taxon.objects.filter(ncbi_taxon_id=int(summaryData["Taxid"])).first()

        level_dict = {v: k for k, v in dict(Assembly.ASSEMBLY_LEVEL).items()}
        type_dict = {v: k for k, v in dict(Assembly.ASSEMBLY_TYPES).items()}

        s = Assembly(type=Resource.RESOURCE_TYPES.ASSEMBLY, name=acc + "_" + name,
                     description=summaryData["AssemblyDescription"],
                     ncbi_tax=tax,
                     ncbi_org=summaryData["SubmitterOrganization"],
                     level=level_dict[summaryData["AssemblyStatus"].lower()],
                     assembly_type=type_dict[summaryData["AssemblyType"].lower()],
                     species_name=summaryData["SpeciesName"])
        s.url = "https://www.ncbi.nlm.nih.gov/assembly/" + acc
        s.type = s.__class__.TYPE
        try:
            s.intraspecific_name = str(
                summaryData["Biosource"]["InfraspeciesList"][0]["Sub_type"]) + " " + \
                                   summaryData["Biosource"]["InfraspeciesList"][0]["Sub_value"]
        except IndexError:
            pass

        try:
            s.release_date = datetime.datetime.strptime(summaryData["SeqReleaseDate"].split(" ")[0],
                                                        "%Y/%m/%d")
        except ValueError:
            pass
        try:
            s.update_date = datetime.datetime.strptime(summaryData["LastUpdateDate"].split(" ")[0],
                                                       "%Y/%m/%d")
        except ValueError:
            pass

        return s

    def save(self, summaryData, ncbi_id):

        s = self.adapt(summaryData)
        s.save(force_insert=True)

        ncbi_org = Organization.objects.get(name="NCBI")
        s.publishers.add(ncbi_org)

        ExternalId(resource=s, organization=ncbi_org,
                   identifier=s.name, type="accession").save(force_insert=True)
        ExternalId(resource=s, organization=ncbi_org,
                   identifier=ncbi_id, type="identifier").save(force_insert=True)
        return s


class NCBIStructureAdapter:

    def fetch(self, ncbi_id):
        return Entrez.read(retry(lambda: Entrez.esummary(db="structure", id=ncbi_id)), validate=False)[0]

    def fetch_list(self, ncbi_ids):
        fetch = Entrez.read(retry(lambda: Entrez.esummary(db="structure", id=ncbi_ids)))
        return fetch

    def adapt(self, summaryData):
        acc = summaryData["PdbAcc"]
        s = Structure(type=Resource.RESOURCE_TYPES.STRUCTURE, name=acc, description=summaryData["PdbDescr"],
                      method=summaryData["ExpMethod"])
        s.type = s.__class__.TYPE

        try:
            s.deposit_date = datetime.datetime.strptime(summaryData["PdbDepositDate"],
                                                        "%Y/%m/%d %H:%M")  # 2017/08/10 00:00
        except ValueError:
            pass
        if "OrganismList" in summaryData and summaryData["OrganismList"]:
            tax = TaxonName.objects.get(name=summaryData["OrganismList"][0])
            if tax:
                tax = tax.taxon
                s.ncbi_tax = tax
        return s

    def save(self, summaryData, ncbi_id):
        ncbi_org = Organization.objects.get(name="NCBI")
        # qs = Structure.objects.filter(external_ids__identifier=acc, external_ids__type="accession")
        # if qs.exists():
        #     return qs.first()

        s = self.adapt(summaryData)
        s.save()
        s.publishers.add(ncbi_org)

        ExternalId(resource=s, organization=ncbi_org,
                   identifier=s.name, type="accession").save()
        ExternalId(resource=s, organization=ncbi_org,
                   identifier=ncbi_id, type="identifier").save()
        return s


class NCBIGDSAdapter:

    def fetch(self, ncbi_id):
        return Entrez.read(retry(lambda: Entrez.esummary(db="gds", id=ncbi_id)))[0]

    def fetch_list(self, ncbi_ids):
        return Entrez.read(retry(lambda: Entrez.esummary(db="gds", id=ncbi_ids)))

    def adapt(self, summaryData) -> Expression:
        acc = str(summaryData["Accession"])

        s = Expression(type=Resource.RESOURCE_TYPES.EXPRESSION, name=acc,
                       description=summaryData["title"] + "." + summaryData["summary"],
                       gdstype=summaryData["gdsType"])
        s.type = s.__class__.TYPE
        if "OrganismList" in summaryData:
            s.ncbi_org = "||".join(summaryData["OrganismList"])
        if "ExpMethod" in summaryData:
            s.method = str(summaryData["ExpMethod"])

        try:
            s.pdat = datetime.datetime.strptime(summaryData["PDAT"], "%Y/%m/%d")  # 2017/08/10 00:00
        except ValueError:
            pass

        tax = TaxonName.objects.get(name=summaryData["taxon"].split(";")[0])
        if tax:
            tax = tax.taxon
            s.ncbi_tax = tax

        return s

    def save(self, summaryData, ncbi_id):

        s = self.adapt(summaryData)
        s.save(force_insert=True)
        ncbi_org = Organization.objects.get(name="NCBI")
        s.publishers.add(ncbi_org)

        ExternalId(resource=s, organization=ncbi_org,
                   identifier=s.name, type="accession").save()
        ExternalId(resource=s, organization=ncbi_org,
                   identifier=ncbi_id, type="identifier").save()
        return s
