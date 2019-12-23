# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from django.urls import reverse


class Organization(models.Model):
    NCBI = "NCBI"
    BOLD = "BOLD"
    SCOPUS = "SCOPUS"
    EBI = "EBI"
    SNDG = "SNDG"
    RENAORG = "RENAORG"

    TYPE = 30

    name = models.CharField(max_length=255)
    url = models.URLField(null=True)
    country = models.CharField(max_length=200, null=True)
    city = models.CharField(max_length=200, null=True)
    scopus_id = models.CharField(max_length=200, null=True)
    scopus_names = models.TextField(null=True)

    deprecated = models.BooleanField(default=False)
    index_updated = models.BooleanField(default=False)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    source = models.ForeignKey("Organization", on_delete=models.SET_NULL, null=True)
    description = models.TextField(null=True)
    main_identifier = models.CharField(max_length=200, null=True)

    class Meta:
        verbose_name_plural = _("Organizations")

    @staticmethod
    def init_orgs():
        if Organization.objects.filter(name=Organization.NCBI).count() == 0:
            Organization.objects.create(name=Organization.NCBI,
                                        description="National Center for Biotechnology Information",
                                        url="https://www.ncbi.nlm.nih.gov/", country="USA")
        if Organization.objects.filter(name=Organization.BOLD).count() == 0:
            Organization.objects.create(name=Organization.BOLD, description="Barcode of Life Data system",
                                        url="http://www.boldsystems.org/", country="International")
        if Organization.objects.filter(name=Organization.SCOPUS).count() == 0:
            Organization.objects.create(name=Organization.SCOPUS, description="Scopus",
                                        url="https://www.scopus.com", country="International")
        if Organization.objects.filter(name=Organization.EBI).count() == 0:
            Organization.objects.create(name=Organization.EBI, description="European Bioinformatics Institute",
                                        url="https://www.ebi.ac.uk/", country="International")
        if Organization.objects.filter(name=Organization.SNDG).count() == 0:
            Organization.objects.create(name=Organization.SNDG, description="Sistema Nacional de Datos Genómicos",
                                        url="http://datos.sndg.mincyt.gob.ar/", country="Argentina")
        if Organization.objects.filter(name=Organization.RENAORG).count() == 0:
            Organization.objects.create(name=Organization.RENAORG,
                                        description="Registro Nacional de las Organizaciones",
                                        url="", country="Argentina")

    def rtype(self):
        return 30  # "org"

    def get_absolute_url(self):
        return reverse('bioresources:organization_view', args=[str(self.id)])

    def targets(self):
        return [rel.external for rel in self.target_rels.all()]

    def __str__(self):
        return " ".join([x for x in [self.name, "|", self.country, self.city] if x])


class OrgRelationship(models.Model):
    TYPES = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            "dependency", "external",
        ])]
    )
    source = models.ForeignKey(Organization, on_delete=models.CASCADE, related_name="target_rels")
    target = models.OneToOneField(Organization, on_delete=models.CASCADE, related_name="source_rels")
    rel_type = models.PositiveIntegerField(choices=TYPES)


class OrgNames(models.Model):
    main = models.ForeignKey(Organization, on_delete=models.CASCADE, related_name="names")
    name = models.CharField(max_length=255)
