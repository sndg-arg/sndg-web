# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from .Resource import Resource


class Assembly(Resource):
    ASSEMBLY_LEVEL = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            "complete genome", 'chromosome', "scaffold", "contig",
        ])]
    )
    ASSEMBLY_TYPES = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            "haploid", 'diploid', "other", "unresolved-diploid", "alternate-pseudohaplotype"
        ])]
    )

    TYPE = Resource.RESOURCE_TYPES.ASSEMBLY

    intraspecific_name = models.CharField(max_length=250, null=True,blank=True)
    species_name = models.CharField(max_length=200, null=True,blank=True)
    level = models.PositiveIntegerField(null=True, choices=ASSEMBLY_LEVEL)

    ncbi_org = models.CharField(max_length=200, null=True,blank=True)
    release_date = models.DateField(null=True)
    update_date = models.DateField(null=True)
    assembly_type = models.PositiveIntegerField(null=True, choices=ASSEMBLY_TYPES)

    class Meta:
        verbose_name_plural = _("Assemblies")

    def related_author_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.author_names()
        return list(set(ran))

    def related_org_names(self):
        ran = []
        for affiliation in self.targets.all():
            ran += affiliation.source.affiliation_names()
        return list(set(ran))
