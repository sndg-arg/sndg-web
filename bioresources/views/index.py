# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as _
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from collections import defaultdict

from haystack.query import SearchQuerySet
from haystack.inputs import  Exact

from ..models.Resource import Resource
from ..models.Tool import Tool
from ..models.Structure import Structure
from ..models.Barcode import Barcode
from ..models.Assembly import Assembly
from ..models.BioProject import BioProject
from ..models.Publication import Publication
from ..models.Expression import Expression
from ..models.ReadsArchive import ReadsArchive


resources = [
    # {"name": "Genomas", "count": 99, "icon": "circle", "type": "genome"},
    # {"name": _("Proteins"), "count": 99,
    #  "icon": "puzzle-piece",
    #  "type": Resource.RESOURCE_TYPES.PROTEIN},
    {"name": _(Tool._meta.verbose_name_plural), "count": 99, "icon": "wrench",
     "type": Resource.RESOURCE_TYPES.TOOL},
    {"name": Structure._meta.verbose_name_plural, "count": 99,
     "icon": "sitemap",
     "type": Resource.RESOURCE_TYPES.STRUCTURE},
    {"name": Barcode._meta.verbose_name_plural, "count": 99,
     "icon": "barcode",
     "type": Resource.RESOURCE_TYPES.BARCODE},
    {"name": Assembly._meta.verbose_name_plural, "count": 99,
     "icon": "adn",
     "type": Resource.RESOURCE_TYPES.ASSEMBLY},
    # {"name": BioProject._meta.verbose_name_plural, "count": 99,
    #  "icon": "briefcase",
    #  "type": Resource.RESOURCE_TYPES.BIOPROJECT},
    {"name": Publication._meta.verbose_name_plural, "count": 99,
     "icon": "book",
     "type": Resource.RESOURCE_TYPES.PUBLICATION},
    {"name": _("Organizations"), "count": 99, "icon": "building", "type": "30"},
    {"name": _("Persons"), "count": 99, "icon": "user", "type": "20"},
    {"name": Expression._meta.verbose_name_plural, "count": 99,
     "icon": "sliders-h",
     "type": Resource.RESOURCE_TYPES.EXPRESSION},
    {"name": ReadsArchive._meta.verbose_name_plural, "count": 99,
     "icon": "sliders-h",
     "type": Resource.RESOURCE_TYPES.READS},

]


def faq(request):
    return render(request, 'faq.html', {})

def index(request):
    db = request.GET.get("db", "all")
    search = request.GET.get("search", "").strip()

    selected = {x: Exact(request.GET[x]) for x in ["authors", "affiliations", "taxon"]
                if x in request.GET}
    if search:
        sqs = SearchQuerySet().filter(content=search, **selected).facet("type")
    else:
        sqs = SearchQuerySet().filter(**selected).facet("type")
        search = "*"

    params = dict(request.GET)
    params["q"] = [search]

    for ft in ["authors", "affiliations", "taxon"]:
        if ft not in selected:
            sqs = sqs.facet(ft, limit=5)

    facets = sqs.facet_counts()
    if "fields" in facets:
        rdata = defaultdict(lambda: 0, {k: v for k, v in facets["fields"]["type"]})
    else:
        rdata = defaultdict(lambda: 0)
    count = 0
    for r in resources:
        r["count"] = rdata[str(r["type"])]
        count += rdata[r["type"]]

    suggestions = []
    if count == 0:
        suggestions = SearchQuerySet().auto_query(search).spelling_suggestion()
        if suggestions:
            suggestions = [x.strip() for x in suggestions.replace("(", " ").split(")") if x.strip()]

    if "fields" in facets:
        del facets["fields"]["type"]
    else:
        facets["fields"] = {}

    return render(request, 'index.html', {
        "stats": resources, "search": search if search != "*" else "", "selected": selected,
        "db": db, "suggestions": suggestions, "querystring": params,
        "sidebarleft": facets["fields"], "sidebarrigth": {"news": [{"title": "n1", "text": "lalala"}]}})

