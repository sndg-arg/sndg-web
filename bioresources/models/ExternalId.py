# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models

from .Organization import Organization
from .Resource import Resource

class ExternalId(models.Model):
    organization = models.ForeignKey(Organization, on_delete=models.CASCADE)
    identifier = models.CharField(max_length=100)
    type = models.CharField(max_length=20)
    resource = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="external_ids")

    class Meta:
        verbose_name_plural = _("ExternalId")

    def __str__(self):
        return self.organization.name + ":" + self.identifier