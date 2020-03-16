# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render
from django.db import transaction
from django.conf import settings
from django.http import HttpResponseNotFound

from bioresources.io.NCBISearch import NCBISearch
from bioresources.models.ExternalId import ExternalId
from bioresources.models.Organization import Organization
from bioresources.models.Structure import Structure
from bioresources.models.Publication import Publication
from bioresources.models.Resource import Collaboration
from bioresources.io.scopus import ScopusDS

from bioresources.io.adapters import scopus_extended_publication
from pdbdb.io.PDB2SQL import PDB2SQL
from bioresources.models.Job import Job
from bioresources.tasks import execute_job
from bioresources.models.jobs.LoadPDBJob import LoadPDBJob
from bioresources.models.jobs.LoadGenomeJob import LoadGenomeFromNCBIJob

from cache_memoize import cache_memoize


class LineResult():
    def __init__(self, line):
        self.line = line
        self.resources = []
        self.message = ""


class ResourceResult():
    def __init__(self, msg, resource, external_id, rtype):
        self.msg = msg
        self.resource = resource
        self.external_id = external_id
        self.rtype = rtype


rtypes = ["assembly", "sra", "gds", "biosample", "structure", "pubmed","bioproject"]

from django.contrib.auth.decorators import login_required
@login_required
def SubmissionImportView(request):
    current_user = request.user
    if request.method == 'POST':

        search = request.POST["search"].strip()
        # rid = request.POST["rid"]
        db = request.POST["rtype"]
        if db == "structure":
            ex = Structure.objects.filter(name__iexact=search)
            if ex.count():
                rr = ResourceResult(__("Already in DB"), ex.get(), search, db)
                results = [rr]
            else:
                results = process_db_id(db, search)

        elif db == "pubmed":
            ex = Publication.objects.filter(doi__icontains=search)
            if ex.count():
                rr = ResourceResult(__("Already in DB"), ex.get(), search, db)
                results = [rr]
            else:
                sds = ScopusDS(settings.SCOPUS_API)
                results = sds.doi(search)
                if results and "error" not in results[0]:
                    doi = results[0]['prism:doi']
                    rr = ResourceResult("", Publication(name=doi, description=results[0]["dc:title"],
                                                        type=Publication.TYPE), doi, db)
                    results = [rr]
                else:
                    results = []
        else:
            if search:
                results = process_db_id(db, search)
            else:
                results = []

        return render(request, 'submission/submission_import.html', {"rtype": db,
                                                                     "collaboration_types": {x[0]: str(x[1]) for x in
                                                                                             Collaboration.COLLABORATION_TYPES},
                                                                     'resource': "assembly", "search": search,
                                                                     "results": results,
                                                                     "count": len(results)})
    else:
        search = request.GET.get("search", "")
        return render(request, 'submission/submission_import.html', {
            'resource': "assembly", "search": search})


def process_search_lines(search):
    results = []
    for l in search.split("\n"):
        if l.strip():
            results.append(process_line(l))
    return results


def process_line(line):
    r = LineResult(line)
    if len(line.split()) > 1:
        r.message = __("lines must not have any spaces")

    if any([x in line for x in
            ["www.ncbi.nlm.nih.gov", "www.ebi.ac.uk", "www.uniprot.org", "rnacentral.org", "ensembl.org"]]):
        return process_url(line)
    elif "http://" in line or "https://" in line:
        r.message = __("url not recognized")

    else:
        r = process_id(line.strip())
    return r


# @cache_memoize(60 * 15)
def process_db_id(db, rid):
    line_result = []
    ncbi_search = NCBISearch()
    r = ExternalId.objects.prefetch_related("resource").filter(resource__type=NCBISearch.db_type[db], type="accession",
                                                               organization=Organization.objects.get(name="NCBI"),
                                                               identifier=rid)
    if r.count():
        r = r.get().resource
        line_result.append(ResourceResult(__("Already in DB"), r,rid,r.type_name()))
    else:
        records, _ = ncbi_search.search_database(db, rid)
        if records:
            for ncbi_id, r in records.items():
                ex = ExternalId.objects.prefetch_related("resource").filter(
                    resource__type=NCBISearch.db_type[db], type="identifier",
                    organization=Organization.objects.get(name="NCBI"), identifier=ncbi_id)
                if ex.count():
                    rr = ResourceResult(__("Already in DB"), ex.get().resource, ncbi_id, db)
                    line_result.append(rr)
                else:
                    rr = ResourceResult("", r, ncbi_id, db)
                    line_result.append(rr)

                #     line_result.message = __("Too many results with this identifier. Is there other bioproject or assembly that groups this resources?")
    return line_result


def process_id(rid):
    line_result = LineResult(rid)
    ncbi_search = NCBISearch()
    counts = ncbi_search.search_all(rid)
    for db, count in counts.items():
        if count:
            if (count == 1) and (db not in ["taxonomy", "genome", "biosample"]):
                process_db_id(db, rid)

    return line_result

from django.contrib.auth.decorators import login_required
@login_required
def SubmitImportView(request):
    ncbi_search = NCBISearch()
    accessions = ["_".join(x.split("_")[1:]) for x in request.POST if x.startswith("accession_")]

    if request.method != 'POST':
        return HttpResponseNotFound()
    for acc in accessions:
        ncbi_id = request.POST["id_" + acc]
        rtype = request.POST["rtype_" + acc]
        relation = int(request.POST["relation_" + acc])

        if rtype == "publication":
            sds = ScopusDS(settings.SCOPUS_API)
            results = sds.doi(ncbi_id)
            record = scopus_extended_publication(results[0])

        else:
            record = ncbi_search.save_resource(rtype, ncbi_id)
            if rtype=="structure":
                with transaction.atomic():
                    job = LoadPDBJob(pdb=record.name)
                    job.save()
                    job.init()
                    job.queue()
                    job.save()
                execute_job.apply_async(args=(job.id,),countdown=10)
            if rtype =="assembly":
                with transaction.atomic():
                    job = LoadGenomeFromNCBIJob(assembly=record)
                    job.save()
                    job.init()
                    job.queue()
                    job.save()
                execute_job.apply_async(args=(job.id,),countdown=10)
        Collaboration(resource=record, person=request.user.person, type=relation).save(force_insert=True)

    return redirect(reverse("bioresources:user_resources"))

# results = process_search_lines(search)
# ncbi_search = NCBISearch()
# page = Page.from_request(request)
# records, total = ncbi_search.search(search, retmax=page.size, retstart=page.offset)
# existing = [x for x in Assembly.objects.filter(name__in=[x.name for x in records])]
#
# if current_user:
#     collaborated = [x.name for x in existing if
#                     hasattr(current_user, "person") and current_user.person in x.collaborators.all()]
# else:
#     collaborated = []
# collaborate = [x.name for x in existing if x.name not in collaborated]
#
# page.set_count(total)
# return render(request, 'submission/import_resource.html',
#               {'resource': "assembly", "search": search, "collaborate": collaborate,
#                "page_obj": page, "results": records, "collaborated": collaborated})
