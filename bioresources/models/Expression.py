# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from .Resource import Resource

class Expression(Resource):
    TYPE = Resource.RESOURCE_TYPES.EXPRESSION

    pdat = models.DateField(null=True,blank=True)
    gdstype = models.CharField(max_length=250, null=True,blank=True,
                               help_text="Tipo de análisis, por ejemplo: 'RNAs no codificantes'")
    submitters = models.TextField(null=True)

    class Meta:
        verbose_name_plural = _("Expression")

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
