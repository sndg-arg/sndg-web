# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models

from .Resource import Resource
from .Organization import Organization
from bioseq.models.Term import Term


class ResourceProperty(models.Model):
    organization = models.ForeignKey(Organization, on_delete=models.DO_NOTHING)
    resource = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="properties")
    term = models.ForeignKey(Term, on_delete=models.DO_NOTHING)

    def __str__(self):
        return str(self.term)

class ResourcePropertyValue(models.Model):
    property = models.OneToOneField(ResourceProperty, on_delete=models.CASCADE, related_name="value")
    value = models.CharField(max_length=200)