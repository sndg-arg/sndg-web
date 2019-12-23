# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioseq.models.Bioentry import Bioentry
from bioseq.models.Seqfeature import Seqfeature
from bioseq.models.Biosequence import Biosequence
from bioresources.models.Assembly import Assembly

from bioseq.io.Pagination import Page


def NucleotideView(request, pk):
    be = Bioentry.objects.prefetch_related("dbxrefs__dbxref", "qualifiers__term")

    be = be.get(bioentry_id=pk)
    #     be = be.get(bioentry_id=pk)
    #     feature_list = be.features.all()
    # else:
    #     be = be.get(bioentry_id=pk)
    qs = be.features.exclude(type_term__identifier__in=["source", "gene"])
    page = Page.from_request(request, qs.count())
    qs = qs.prefetch_related(
        "dbxrefs__dbxref", "qualifiers__term", "locations", "source_term", "type_term").all()
    feature_list = [f for f in qs[page.offset():page.end()]
                    ]  # if f.type_term.identifier not in ["source", "gene"]

    assembly = Assembly.objects.get(name=be.biodatabase.name)

    return render(request, 'resources/nucleotide_view.html', {"query": "", "page_obj": page,
                                                              "object": be, "assembly": assembly,
                                                              "feature_list": feature_list,
                                                              "sidebarleft": 0})  # "sidebarrigth": {"news": [{"title": "n1", "text": "lalala"}]
