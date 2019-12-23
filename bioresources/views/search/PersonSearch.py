# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from haystack.generic_views import SearchView
from haystack.inputs import Exact
from haystack.query import SearchQuerySet


def person(request, pk):
    person = Person.objects.get(id=pk)

    sqs = SearchQuerySet().filter(authors=Exact(person.complete_name())).facet("type")

    facets = sqs.facet_counts()
    facets = {k: v for k, v in facets["fields"]["type"] if v > 0}

    search_url = (reverse('bioresources:search_view') +
                  "?q=*&authors=" + person.complete_name() + "&type=")

    return render(request, 'resources/person.html', {
        "person": person, "facets": facets, "search_url": search_url,
        "nameMap": {"person": Person._meta.verbose_name_plural,
                    Resource.RESOURCE_TYPES.PUBLICATION: Publication._meta.verbose_name_plural},
        "sidebarleft": 1, })