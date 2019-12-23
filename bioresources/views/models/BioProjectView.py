# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.BioProject import BioProject
from bioresources.views import labelize

from bioresources.io.GraphRepo import GraphRepo
from bioresources.io.NCBISearch import NCBISearch


def bioproject(request, pk):
    bioproject = BioProject.objects.get(id=pk)

    graph, related_resources = GraphRepo.get_neighborhood(bioproject.id, "BioProject", 1)

    external_ids = [x.identifier for x in bioproject.external_ids.all() if x.type == "accession"]
    external_url = ""
    if external_ids:
        external_url = ("https://www.ncbi.nlm.nih.gov/" + NCBISearch.rtype2ncbb[BioProject.TYPE] + "/" + external_ids[
            0])
    collaboration = request.user.get_collaboration(bioproject) if request.user.is_authenticated else None
    params = {"external_url": external_url,
              "bioproject": bioproject, "graph": graph,
              "related_resources": related_resources,
              "sidebarleft": 1, "collaboration": collaboration}
    return render(request, 'resources/bioproject.html', params)
