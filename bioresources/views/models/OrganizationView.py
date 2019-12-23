# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.Organization import Organization
from bioresources.io.GraphRepo import GraphRepo


def organization(request, pk):
    org = Organization.objects.get(id=pk)
    graph, related_resources = GraphRepo.get_neighborhood(pk, "Organization",level=2)

    return render(request, 'resources/organization.html',
                  {"graph": graph, "related_resources": related_resources, "pk": pk,
                   "rtype_src":"Organization","organization": org, "sidebarleft": 1, })
