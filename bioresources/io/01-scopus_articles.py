import datetime
import json
import os

from tqdm import tqdm

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "sndg.settings")
import django

django.setup()

from bioresources.models import Publication, Organization, Person, Resource, Affiliation

with open("/home/eze/workspace/genomica-arg/data/scopus3.json") as h:
    articles = json.load(h)

# publications = []
err1 = 0
err2 = 0
pepe = []
for i,article in enumerate(tqdm(articles)):

    if Publication.objects.filter(scopus_id=article['dc:identifier']).exists():
        publication = Publication.objects.get(scopus_id=article['dc:identifier'])
    elif ("prism:doi" in article) and \
            Publication.objects.filter(doi=article['prism:doi']).exists():
        publication = Publication.objects.get(doi=article['prism:doi'])
    elif Publication.objects.filter(name=article["dc:title"][:350]).exists():
        publication = Publication.objects.get(name=article["dc:title"][:350])
    else:

        publication = Publication(
            type=Resource.PUBLICATION,
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
        publication.save()
    arg = []
    for affiliation in article["affiliation"]:

        afcountry = affiliation["affiliation-country"]

        org = Organization(name=affiliation["affilname"],
                           country=afcountry,
                           city=affiliation["affiliation-city"]
                           )
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

        if org.country == "Argentina" and org.scopus_id:
            arg.append(org.scopus_id)

    if not arg:
        err2 += 1
        for a in article["affiliation"]:
            if not a["affiliation-country"]:
                if a["affiliation-city"]:
                    pepe.append(a["affiliation-city"])
        continue

    if "author" in article:
        for author in article["author"]:
            person = Person(surname=author["surname"],
                            name=author["given-name"] if author["given-name"] else "",
                            scopus_id=author["authid"])
            if Person.objects.filter(scopus_id=person.scopus_id).exists():
                person = Person.objects.get(scopus_id=person.scopus_id)
            else:
                person.save()

            if ("afid" in author):

                aff = Affiliation(publication=publication, author=person)
                aff.save()
                for affdict in author["afid"]:
                    aff.organizations.add(Organization.objects.get(scopus_id=affdict["$"]))
                aff.save()

                if [x for x in author["afid"] if x["$"] in arg]:
                    # person.arg_affiliation = True
                    person.save()

print([err2])
print(sorted(list(set(pepe))))
# Publication.objects.bulk_create(publications)
