# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models

from .Resource import Resource

class Sample(Resource):
    TYPE = Resource.RESOURCE_TYPES.SAMPLE

    origin_props = ['Origin (developed or donated from)', 'geo_loc_name', 'geographic location (country and/or sea)',
                    'birth_location', 'country_of_birth', 'geo-loc-name']

    country = models.CharField(max_length=100)
    subdivision = models.CharField(max_length=150)
    collection_date = models.DateField(null=True)
    publication_date = models.DateField(null=True)
    update_date = models.DateField(null=True)

    lat = models.FloatField(null=True,blank=True)
    lon = models.FloatField(null=True,blank=True)

    class Meta:
        verbose_name_plural = _("Samples")

    def origin_dict(self):
        props = self.properties.prefetch_related("value").filter(
            Q(term__ontology__name="NCBI sample") & Q(term__identifier__in=Sample.origin_props))
        return {prop.term.name: prop.value.value for prop in props}