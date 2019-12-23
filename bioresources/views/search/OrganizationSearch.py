# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from haystack.generic_views import SearchView
from haystack.inputs import Exact
from haystack.query import SearchQuerySet


def organization(request, pk):
    organization = Organization.objects.get(id=pk)

    sqs = SearchQuerySet().filter(affiliations=Exact(organization.name)).facet("type")

    facets = sqs.facet_counts()
    facets = {k: v for k, v in facets["fields"]["type"] if v > 0}

    search_url = (reverse('bioresources:search_view') +
                  "?q=*&affiliations=" + organization.name + "&type=")
    return render(request, 'resources/organization.html', {
        "organization": organization, "facets": facets, "search_url": search_url,
        "sidebarleft": 1, "nameMap":
            {"person": Person._meta.verbose_name_plural,
             Resource.RESOURCE_TYPES.PUBLICATION: Publication._meta.verbose_name_plural}})