# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __

from haystack.generic_views import SearchView
from haystack.inputs import Exact
from haystack.query import SearchQuerySet

from django import forms
from haystack.forms import SearchForm
from haystack.generic_views import SearchView

from bioseq.io.Pagination import Page
from bioresources.models.Resource import Resource


class BioSearchForm(SearchForm):
    models = [Resource]

    type = forms.TextInput()
    authors = forms.TextInput()
    affiliations = forms.TextInput()
    taxon = forms.TextInput()

    #
    # # start_date = forms.DateField(required=False)
    # # end_date = forms.DateField(required=False)

    def search(self):
        # First, store the SearchQuerySet received from other processing.
        sqs = super(BioSearchForm, self)

        if not self.is_valid():
            return self.no_query_found()

        if not self.cleaned_data.get('q'):
            return self.no_query_found()
        from haystack.inputs import AutoQuery, Exact

        params = {"content": AutoQuery(self.cleaned_data['q']),
                  "type": Exact(self.data["type"])}
        for x in ["authors", "affiliations", "taxon"] + Resource.facet_dict.get(self.data["type"], []):
            if x in self.data and self.data[x]:
                params[x] = self.data[x]
        sqs = self.searchqueryset.filter(**params)

        if self.load_all:
            sqs = sqs.load_all()

        # # Check to see if a start_date was chosen.
        # if self.cleaned_data['start_date']:
        #     sqs = sqs.filter(pub_date__gte=self.cleaned_data['start_date'])
        #
        # # Check to see if an end_date was chosen.
        # if self.cleaned_data['end_date']:
        #     sqs = sqs.filter(pub_date__lte=self.cleaned_data['end_date'])

        return sqs

class BioSearchView(SearchView):


    form_class = BioSearchForm

    def get_queryset(self):
        queryset = super(BioSearchView, self).get_queryset()
        rtype = self.request.GET["type"]

        for x in ["authors", "affiliations", "taxon"] + Resource.facet_dict.get(rtype, []):
            if x not in self.request.GET:
                queryset = queryset.facet(x, limit=5)

        return queryset  # .filter(pub_date__gte=date(2015, 1, 1))

    def form_valid(self,form):
        if not form.cleaned_data["q"].strip():
            form.cleaned_data["q"] = "*"

        return super(BioSearchView, self).form_valid(form)


    def get_context_data(self, *args, **kwargs):
        if kwargs["query"].strip() == "*":
            kwargs["query"] = ""
        context = super(BioSearchView, self).get_context_data(*args, **kwargs)

        context["sidebarleft"] = {k: [(y1, y2) for y1, y2 in v if y2] for k, v in
                                  context["view"].queryset.query._facet_counts["fields"].items()}

        if hasattr(context["view"].queryset.query, "_raw"):
            context["suggestions"] = context["view"].queryset.query._raw.spellcheck["suggestions"]
        # dbtitle = self.request.GET["db"]
        # type
        # authors
        # affiliations
        context["dbtype"] = self.get_form_kwargs()["data"]["type"]

        prange = list(context["page_obj"].paginator.page_range)

        context["page_obj"].paginator.show_pages = (prange[max(
            context["page_obj"].number - 3, 0):context["page_obj"].number + 2])

        context["selected"] = {x: self.get_form_kwargs()["data"][x] for x in
                               ["authors", "affiliations", "taxon"] +
                               Resource.facet_dict.get(self.request.GET["type"], [])
                               if x in self.get_form_kwargs()["data"]}

        context["params"] = dict(self.request.GET)
        suggestion_list = []
        context["suggestion_list"] = suggestion_list
        if "suggestions" in context and context["suggestions"]:
            for x in context["suggestions"][1::2]:
                for y in x["suggestion"]:
                    suggestion_list.append(y)

        # SearchQuerySet().filter(content="genomea") .facet('affiliations',limit=5).facet("type").facet_counts()
        return context
