# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from .Resource import Resource

class ReadsArchive(Resource):
    TYPE = Resource.RESOURCE_TYPES.READS

    release_date = models.DateField(null=True)
    update_date = models.DateField(null=True)

    class Meta:
        verbose_name_plural = _("Reads Archive")