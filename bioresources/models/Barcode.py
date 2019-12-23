# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from .Resource import Resource

class Barcode(Resource):
    TYPE = Resource.RESOURCE_TYPES.BARCODE
    BIODBNAME = "Barcodes"

    country = models.CharField(max_length=100)
    subdivision = models.CharField(max_length=150)
    marker = models.CharField(max_length=50, null=True)
    image_url = models.URLField(null=True)
    bold_org = models.CharField(max_length=255, null=True)
    collectors = models.CharField(max_length=255, null=True)

    lat = models.FloatField(null=True,blank=True)
    lon = models.FloatField(null=True,blank=True)

    class Meta:
        verbose_name_plural = _("Barcodes")