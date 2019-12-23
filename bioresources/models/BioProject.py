# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from .Resource import Resource

class BioProject(Resource):
    """
    https://www.ncbi.nlm.nih.gov/books/NBK169438/
    """

    SAMPLE_SCOPE_TYPES = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            "monoisolate", 'multi-species', "environment", "synthetic", "other",
        ])]
    )

    MATERIAL_TYPES = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            "genome", 'metagenome', "chromosome", "transcriptome", "reagent", "proteome",
        ])])

    CAPTURE_TYPES = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            "whole", 'exome', "barcode", "TargetedLocusLoci",
        ])]
    )

    TYPE = Resource.RESOURCE_TYPES.BIOPROJECT

    sample_scope = models.PositiveIntegerField(choices=SAMPLE_SCOPE_TYPES, null=True,blank=True)
    material = models.PositiveIntegerField(choices=MATERIAL_TYPES, null=True,blank=True)
    capture = models.PositiveIntegerField(choices=CAPTURE_TYPES, null=True,blank=True)

    target = models.CharField(max_length=200, null=True)
    submitters = models.TextField(null=True)

    method = models.CharField(max_length=200, null=True)

    # objetive = models.CharField(max_length=200, blank=False)

    class Meta:
        verbose_name_plural = _("BioProjects")

    def related_author_names(self):
        ran = []
        for affiliation in self.targets.filter(source__type=Resource.RESOURCE_TYPES.PUBLICATION).all():
            ran += affiliation.source.author_names()
        return list(set(ran))

    def related_org_names(self):
        ran = []
        for affiliation in self.targets.filter(source__type=Resource.RESOURCE_TYPES.PUBLICATION).all():
            ran += affiliation.source.affiliation_names()
        return list(set(ran))
