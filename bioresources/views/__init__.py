# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from collections import defaultdict

from bioseq.io.Pagination import Page

from bioseq.views import labelize

# from django.forms import Form
from django.http import HttpResponseRedirect, HttpResponse

# from .filters import PublicationFilter
# from .forms import BioSearchForm, AssemblyForm


# from django_filters.views import FilterView
# from django_tables2.views import SingleTableMixin


# class FilteredPersonListView(SingleTableMixin, FilterView):
#     table_class = PublicationTable
#     model = Publication
#     template_name = 'publications/list.html'
#
#     filterset_class = PublicationFilter
#     paginate = {'per_page': 25}


# def publications(request):
#     table = PublicationTable(Publication.objects.all())
#     RequestConfig(request,paginate={'per_page': 25}).configure(table)
#     return render(request, 'publications/list.html',
#                   {"publications": table})

from bioseq.models.Taxon import TaxIdx


def tax_data(request):
    qs = TaxIdx.objects.filter(text__icontains=request.GET["q"])
    tot = qs.count()

    return {
        "results": [{"id": x.tax_id, "text": x.tax.scientific_name()} for x in qs[:10]]
    }
