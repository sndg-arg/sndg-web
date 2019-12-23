# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.ReadsArchive import ReadsArchive
from bioresources.io.NCBISearch import NCBISearch
from bioresources.io.GraphRepo import GraphRepo


def reads(request, pk):
    sra = ReadsArchive.objects.prefetch_related("external_ids").get(id=pk)

    graph, related_resources = GraphRepo.get_neighborhood(pk, "Reads", level=1)

    external_ids = [x.identifier for x in sra.external_ids.all() if x.type == "accession"]

    external_url = ""
    if external_ids:
        external_url = ("https://www.ncbi.nlm.nih.gov/" + NCBISearch.rtype2ncbb[ReadsArchive.TYPE] + "/" + external_ids[
            0] + "[accn]")

    collaboration = request.user.get_collaboration(sra) if request.user.is_authenticated else None
    params = {"readsarchive": sra, "sidebarleft": 1, "level": 1, "external_url": external_url,
              "graph": graph, "related_resources": related_resources, "pk": pk, "rtype_src": "Reads",
              "collaboration": collaboration}
    return render(request, 'resources/readsarchive.html', params)
