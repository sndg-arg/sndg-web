# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.io.GraphRepo import GraphRepo
from bioresources.models.Sample import Sample


def sample_view(request, pk):
    sample = Sample.objects.get(id=pk)
    graph, related_resources = GraphRepo.get_neighborhood(pk, "Sample", level=2)

    return render(request, 'resources/sample.html', {
        "graph": graph, "related_resources": related_resources, "pk": pk, "rtype_src": "Sample",
        "sample": sample, "sidebarleft": 1,"level":1})
