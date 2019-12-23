from django.db import models
from django.db.models import Prefetch, Q
from django.db.models.query import QuerySet

from polymorphic.managers import PolymorphicManager
from polymorphic.managers import PolymorphicQuerySet
from polymorphic.managers import PolymorphicManager

class ResourceQuerySet(PolymorphicQuerySet):

    def oai_compliant(self):
        # return self.exclude(Q(creators=None) | Q(publishers=None)).prefetch_related(
        #     "creators", "publishers"
        # ).filter(deprecated=False, index_updated=False)
        return self.filter(deprecated=False, index_updated=False)

    def publication_related(self, country, only_not_indexed=False):
        from ..models.Publication import Publication
        qs = (Publication.objects.prefetch_related(
            "affiliations__organizations", "affiliations__author")
            .filter(affiliations__organizations__country=country))
        q = {
            "deprecated": False
        }
        if only_not_indexed:
            q["index_updated"] = False
        return (self.prefetch_related("ncbi_tax__names")
            .prefetch_related(Prefetch("targets__source",
                                       queryset=qs))
            .filter(**q))


class BioResourceManager(PolymorphicManager):

    def get_query_set(self):
        return ResourceQuerySet(self.model)

    def __getattr__(self, attr, *args):
        # see https://code.djangoproject.com/ticket/15062 for details
        if attr.startswith("_"):
            raise AttributeError
        return getattr(self.get_query_set(), attr, *args)
