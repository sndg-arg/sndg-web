# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.Structure import Structure
from bioresources.io.GraphRepo import GraphRepo
from bioresources.models.jobs.LoadPDBJob import LoadPDBJob
from bioresources.models.Job import Job

def structure(request, pk):
    pdb = Structure.objects.get(id=pk)
    graph, related_resources = GraphRepo.get_neighborhood(pk, "Structure", level=1)
    external_url = "https://www.rcsb.org/structure/" + pdb.name

    loaded = True
    job = LoadPDBJob.objects.filter(pdb=pdb.name)
    if job.exists():
        job = job.order_by("-id")[0]
        if job.status != Job.STATUS.FINISHED:
            loaded = False

    collaboration = request.user.get_collaboration(pdb) if request.user.is_authenticated else None
    params =  {"external_url":external_url,"loaded":loaded,"collaboration":collaboration,
               "graph": graph, "related_resources": related_resources, "pk": pk, "rtype_src": "Structure",
               "pdb": pdb, "sidebarleft": 1, "level": 1}
    return render(request, 'resources/structure.html',params)
