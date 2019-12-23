# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models

from .Resource import Resource


class ResourceRelation(models.Model):
    source = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="targets")
    target = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="sources")
    role = models.CharField(max_length=200, blank=False)
    deprecated = models.BooleanField(default=False)

    class Meta:
        unique_together = (('source', 'target', 'deprecated'),)
        verbose_name_plural = _("Resource Relations")

    def __str__(self):
        return ("(" + str(self.source.type) + ":" + str(self.source.id) +
                ") -> " + "(" + str(self.target.type) + ":" + str(self.target.id) + ")")