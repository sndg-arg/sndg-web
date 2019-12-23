import json
import os
from datetime import datetime

from tqdm import tqdm

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "sndg.settings")
import django

django.setup()

from bioresources.models import Assembly, ExternalId, Organization

with open("/home/eze/workspace/genomica-arg/data/scopus_ncbi_assemblies.json2") as h:
    lineas = [x for x in h.readlines() if x.strip()]

ncbi = Organization.objects.get(name="NCBI")

for l in tqdm(lineas):
    data = json.loads(l)
    if not Assembly.objects.filter(name=data["acc"] + " " + data["name"]).exists():
        s = Assembly(type=data["db"], name=data["acc"] + " " + data["name"],
                     description=data["acc"] + " " + data["name"],

                     species_name=data["SpeciesName"],
                     level=data["AssemblyStatus"],
                     ncbi_org=data["org"],
                     assembly_type=data["AssemblyType"] if data["AssemblyType"] in ["haploid", "diploid"] else "other")
        s.ncbi_tax = int(data["taxid"])
        if "intraspecific_name" in data:
            s.intraspecific_name = data["intraspecific_name"]

        try:
            s.release_date = datetime.strptime(data["release"].split(" ")[0], "%Y/%m/%d")  # 2017/08/10 00:00
        except:
            pass
        try:
            s.update_date = datetime.strptime(data["update"].split(" ")[0], "%Y/%m/%d")  # 2017/08/10 00:00
        except:
            pass

        s.save()

        eid = ExternalId(organization=ncbi, identifier=data["assembly_id"],
                         type="identifier", resource=s)
        eid.save()
        eid = ExternalId(organization=ncbi, identifier=data["acc"],
                         type="accession", resource=s)
        eid.save()
