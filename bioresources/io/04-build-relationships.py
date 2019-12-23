import json
import os

from tqdm import tqdm

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "sndg.settings")
import django

django.setup()

from bioresources.models import BioProject, Assembly, Structure, Expression, ResourceRelation

mappings = {}
assembly_nuc = {}
dbmap = {"gds": Expression, "bioproject": BioProject, "nuccore": Assembly, "structure": Structure}

with open("/home/eze/workspace/genomica-arg/data/scopus_ncbi_assemblies.json2") as h:
    lineas = [x for x in h.readlines() if x.strip()]
    for l in tqdm(lineas):
        data = json.loads(l)
        if data["from_db"] == "nuccore":
            assembly_nuc[data["from_id"]] = data["assembly_id"]
        else:
            t = Assembly.objects.get(external_ids__identifier=data["acc"])
            dbmap = {"gds": Expression, "bioproject": BioProject, "nuccore": Assembly, "structure": Structure}
            for p in dbmap[data["from_db"]].objects.filter(external_ids__identifier=data["from_id"],
                                                external_ids__type="identifier"):
                ResourceRelation.objects.get_or_create(source =p, target =t, role=data["from_db"] + "_assembly" )


# df = pd.read_table("/home/eze/workspace/genomica-arg/data/scopus_ncbi_links.txt", sep="\t",
#                    names=["doi", "pmc", "ncbi_id", "link", "ids"], index_col=False)
#
# for _, r in df.iterrows():
#     k = r["link"].split("_")[1] + "|" + str(r["doi"] + "|" + str(r["ncbi_id"]))
#     if r["link"].split("_")[1] == "nuccore":
#         seqid = r["ids"].split(",")[0]
#         if seqid in assembly_nuc:
#             mappings[k] = assembly_nuc[seqid]
#         else:
#             print("seq sin assembly...")
#
#     else:
#         mappings[k] = r["ids"].split(",")[0]
#
#
# for k, v in tqdm(mappings.items()):
#
#     db, doi, ncbi_id = k.split("|")
#     if db == "sra":
#         continue
#     try:
#         p = Publication.objects.get(doi=doi)
#     except ObjectDoesNotExist:
#         try:
#             p = Publication.objects.get(pubmed_id=ncbi_id)
#         except:
#             print ((doi,ncbi_id) )
#
#
#     rs = list(dbmap[db].objects.filter(external_ids__identifier=v, external_ids__type="identifier"))
#     if not rs:
#         print ("NO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1")
#     for r in rs:
#         ResourceRelation.objects.get_or_create(source =p, target =r, role=p.type + "_" + r.type)
