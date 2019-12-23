# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.Publication import Publication
from bioresources.io.GraphRepo import GraphRepo

from bioresources.models.Affiliation import Affiliation
from bioresources.models.Organization import Organization

def publication(request, pk):




    publication = Publication.objects.prefetch_related("targets__target", "affiliations__author",
                                                       "affiliations__organizations").get(id=pk)
    orgs = publication.affiliation_names()
    org_map = {}
    authors = []
    for aff in publication.affiliations.all():
        author = aff.author
        author.affs = " ".join(
            ["(" + str(orgs.index(org.name) + 1) + ")" for org in aff.organizations.all()])
        if author not in authors:
            authors.append(author)
        for org in aff.organizations.all():
            org_map[org.name] = org

    graph, related_resources = GraphRepo.get_neighborhood(pk, "Publication", level=1)

    return render(request, 'resources/publication.html', {
        "graph": graph, "related_resources": related_resources, "pk": pk, "rtype_src": "Publication", "level": 1,
        "publication": publication, "authors": authors, "org_map": org_map, "sidebarleft": 1})
