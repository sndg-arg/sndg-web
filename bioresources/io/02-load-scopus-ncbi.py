import json
import os
from datetime import datetime

from tqdm import tqdm

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "sndg.settings")
import django

django.setup()

from bioresources.models import BioProject, Structure, Expression, ExternalId, Organization
from biosql.models import TaxonName

with open("/home/eze/workspace/genomica-arg/data/scopus_ncbi_data2.json") as h:
    lineas = [x for x in h.readlines() if x.strip()]

ncbi = Organization.objects.get(name="NCBI")

for l in tqdm(lineas):
    data = json.loads(l)

    if data["db"] == "structure":
        if not Structure.objects.filter(name=data["name"]).exists():
            s = Structure(type=data["db"], name=data["name"], description=data["description"],
                          pdbClass=data["pdbClass"], method=data["method"], org_list=data["org_list"])
            try:
                s.deposit_date = datetime.strptime(data["deposit_date"], "%Y/%m/%d %H:%M")  # 2017/08/10 00:00
            except:
                pass

            tax = TaxonName.objects.get(name=data["org_list"].split("|")[0])
            if tax:
                tax = tax.taxon.ncbi_taxon_id
                s.ncbi_tax = tax

            s.save()

            eid = ExternalId(organization=ncbi, identifier=data["identifier"],
                             type="identifier", resource=s)
            eid.save()

    elif data["db"] == "gds":
        if not Expression.objects.filter(name=data["name"]).exists():
            s = Expression(type=data["db"], name=data["name"], description=data["description"],
                           gdstype=data["gdstype"], submitters=data["submitters"])
            try:
                s.pdat = datetime.strptime(data["pdat"], "%Y/%m/%d")  # 2017/08/10 00:00
            except:
                pass

            tax = TaxonName.objects.get(name=data["tax"].split(";")[0])
            if tax:
                tax = tax.taxon.ncbi_taxon_id
                s.ncbi_tax = tax

            s.save()

            eid = ExternalId(organization=ncbi, identifier=data["identifier"],
                             type="identifier", resource=s)
            eid.save()
            eid = ExternalId(organization=ncbi, identifier=data["acc"],
                             type="accession", resource=s)
            eid.save()



    elif data["db"] == "bioproject":
        if not BioProject.objects.filter(name=data["name"]).exists():
            s = BioProject(type=data["db"], name=data["name"], description=data["description"],
                           sample_scope=data["scope"], target=data["target"], submitters=data["submitters"],
                           ncbi_tax=int(data ["taxid"]))



            s.save()

            eid = ExternalId(organization=ncbi, identifier=data["identifier"],
                             type="identifier", resource=s)
            eid.save()
            eid = ExternalId(organization=ncbi, identifier=data["acc"],
                             type="accession", resource=s)
            eid.save()









