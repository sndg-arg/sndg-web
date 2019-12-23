from collections import defaultdict
import json
import requests
from datetime import datetime

from bioresources.models.Person import Person
from bioresources.models.Organization import Organization
from bioresources.models.Affiliation import Affiliation
from bioresources.models.Publication import Publication
from bioresources.models.Resource import Collaboration


class EBISearch():
    url_base = "https://www.ebi.ac.uk/europepmc/webservices/rest/"

    @staticmethod
    def doi(doi):
        ebi = Organization.objects.get(name=Organization.EBI)
        url = EBISearch.url_base + "search?resultType=core&format=json&query=doi:" + doi
        res = requests.get(url)
        if res.status_code == 200:
            data = res.json()
            if data["hitCount"]:
                data = data["resultList"]["result"][0]
                persons = []
                affs = defaultdict(lambda: [])
                for i, author in enumerate(data["authorList"]["author"]):
                    p = Person(name=author["firstName"], surname=author["lastName"], source=ebi)
                    persons.append(p)
                    if "affiliation" in author:
                        affs[author["affiliation"].strip()].append(i)
                orgs_list = [data["affiliation"]] if ("affiliation" in data) else []
                orgs = [Organization(name=x, source=ebi) for x in list(affs.keys()) + orgs_list ]

                record = {"doi": data["doi"], "title": data["title"], "abstract": data["abstractText"] if "abstractText" in data else "",
                          "persons": persons, "orgs": orgs, "affs": dict(affs), "record": data}
                return record
            else:
                return None
        else:
            raise Exception("can't connecto to EBI")

    @staticmethod
    def save(record):
        ebi = Organization.objects.get(name=Organization.EBI)
        affs = []
        publication = Publication(name=record["title"], description=record["abstract"], doi=record["doi"],
                                  date_of_publication=datetime.strptime(record["record"]["firstPublicationDate"],
                                                                        "%Y-%m-%d"), )
        publication.save()
        for p in record["persons"]:
            data = {x: p.__dict__[x] for x in ["name", "surname"]}
            data["source"] = ebi
            person = Person.objects.get_or_create(**data)[0]

            aff = Affiliation.objects.get_or_create(
                resource=publication, author=person)[0]
            affs.append(aff)

        for o in record["orgs"]:
            data = {"name": o.name, "source": ebi, "country": 'Argentina' if 'argentina' in  o.name.lower() else '' }

            o = Organization.objects.get_or_create(**data)[0]
            if o.name in record["affs"]:
                for idx in record["affs"][o.name]:
                    affs[idx].organizations.add(o)
            else:
                if not Collaboration.objects.filter(resource=publication,organization=o,
                                                type=Collaboration.COLLABORATION_TYPES.other,info="affiliated").exists():
                    Collaboration.objects.create(resource=publication,organization=o,
                                                 type=Collaboration.COLLABORATION_TYPES.other,info="affiliated")
        return publication
