# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.Tool import Tool
from bioresources.models.Resource import Collaboration,Resource

from bioresources.io.GraphRepo import GraphRepo


def tool(request, pk):
    relate_with = None
    if "relate_with" in request.session:
        relate_with = Resource.objects.get(id=int(request.session["relate_with"]))

    tool = Tool.objects.get(id=pk)
    graph, related_resources = GraphRepo.get_neighborhood(pk, "Tool", level=1)
    collaboration = request.user.get_collaboration(tool)   if request.user.is_authenticated else None

    params = {"external_url": tool.url, "collaboration": collaboration,
              "graph": graph, "related_resources": related_resources, "pk": pk,
              "rtype_src": "Tool",
              "tool": tool, "sidebarleft": 1, "relate_with": relate_with}
    return render(request, 'resources/tool.html', params)
