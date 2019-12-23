# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioseq.views import labelize

from bioresources.models.Expression import Expression
from bioresources.io.GraphRepo import GraphRepo
from bioresources.io.NCBISearch import NCBISearch


def expression(request, pk):
    expression = Expression.objects.get(id=pk)
    graph, related_resources = GraphRepo.get_neighborhood(pk, "Expression", level=2)

    external_ids = [x.identifier for x in expression.external_ids.all() if x.type == "accession"]
    external_url = ""
    if external_ids:
        external_url = ("https://www.ncbi.nlm.nih.gov/" + NCBISearch.rtype2ncbb[Expression.TYPE] + "/" + external_ids[
            0])

    collaboration = request.user.get_collaboration(expression) if request.user.is_authenticated else None
    params = {
        "expression": expression, "sidebarleft": 1, "external_url": external_url, "collaboration": collaboration,
        "graph": graph, "related_resources": related_resources, "pk": pk, "rtype_src": "Expression"
    }

    return render(request, 'resources/expression.html', params)
