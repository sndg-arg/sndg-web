# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from .Resource import Resource

class Publication(Resource):
    TYPE = Resource.RESOURCE_TYPES.PUBLICATION

    doi = models.CharField(max_length=100)
    date_of_publication = models.DateField(max_length=200)
    pubmed_id = models.CharField(max_length=50, null=True)
    electronic_id = models.CharField(max_length=50, null=True)
    scopus_id = models.CharField(max_length=50, null=True)
    issn = models.CharField(max_length=50, null=True)

    associated = models.BooleanField(default=False)

    class Meta:
        verbose_name_plural = _("Publications")

    def affiliation_names(self, country=False):
        affs = []
        for affiliation in self.affiliations.all():
            qs = affiliation.organizations
            if country:
                qs = qs.filter(country=country)
            for org in qs.all():
                if org.name not in affs:
                    affs.append(org.name)
        return affs

    def author_names(self, country=None):
        authors = []
        if country:
            sqs = self.affiliations.filter(organizations__country=country)
        else:
            sqs = self.affiliations

        for affiliation in sqs.all():
            if affiliation.author.complete_name() not in authors:
                authors.append(affiliation.author.complete_name())
        return authors

    def author_names_idx(self, country=None):
        return " ".join(self.author_names(country))

    def affiliation_names_idx(self, country=None):
        return " ".join(self.affiliation_names(country))