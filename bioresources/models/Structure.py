# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models

from .Resource import Resource

class Structure(Resource):
    TYPE = Resource.RESOURCE_TYPES.STRUCTURE

    pdbClass = models.CharField(max_length=50, null=True)
    deposit_date = models.DateField(null=True)
    method = models.CharField(max_length=50, null=True)
    org_list = models.TextField(null=True)

    class Meta:
        verbose_name_plural = _("Structures")

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