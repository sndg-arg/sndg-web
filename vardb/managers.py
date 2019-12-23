"""

"""
from django.db import models
from django.db.models import IntegerField
from django.db.models.functions import Cast
from django.db.models.query import QuerySet


class VariantannotationQuerySet(QuerySet):

    def dp_filter(self, variant_collection, depth=30):
        from .models import Variantannotation
        return Variantannotation.objects.select_related("assignment_fk").filter(
            assignment_fk__variant_collection_fk=variant_collection, prop=Variantannotation.DP).annotate(
            value_int=Cast('value', IntegerField())).filter(value_int__lt=depth)


class VariantannotationManager(models.Manager):

    def get_query_set(self):
        return VariantannotationQuerySet(self.model)

    def __getattr__(self, attr, *args):
        # see https://code.djangoproject.com/ticket/15062 for details
        if attr.startswith("_"):
            raise AttributeError
        return getattr(self.get_query_set(), attr, *args)


class ReportedAlleleQuerySet(QuerySet):

    def exact_reported(self, phenotype, variant_collection):
        from .models import ReportedAllele
        return ReportedAllele.objects.filter(phenotype=phenotype,
                                             effect__alleles__allele_fk__assignments__variant_collection_fk=variant_collection)

    def pos_reported(self, phenotype, variant_collection):
        from .models import ReportedAllele
        return ReportedAllele.objects.filter(phenotype=phenotype,
                                             allele__variant_fk__alleles__assignments__variant_collection_fk=variant_collection)


class ReportedAlleleManager(models.Manager):

    def get_query_set(self):
        return ReportedAlleleQuerySet(self.model)

    def __getattr__(self, attr, *args):
        # see https://code.djangoproject.com/ticket/15062 for details
        if attr.startswith("_"):
            raise AttributeError
        return getattr(self.get_query_set(), attr, *args)


class VariantassignmentQuerySet(QuerySet):

    def reported_positions(self, variant_collection, gene_list):
        from .models import Variantassignment
        return Variantassignment.objects.filter(
            variant_collection_fk=variant_collection,
            variant_fk__gene__in=gene_list)


class VariantassignmentManager(models.Manager):

    def get_query_set(self):
        return VariantassignmentQuerySet(self.model)

    def __getattr__(self, attr, *args):
        # see https://code.djangoproject.com/ticket/15062 for details
        if attr.startswith("_"):
            raise AttributeError
        return getattr(self.get_query_set(), attr, *args)
