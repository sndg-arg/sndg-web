# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.db import models
from django.shortcuts import reverse

class Biodatabase(models.Model):
    biodatabase_id = models.AutoField(primary_key=True)
    name = models.CharField(unique=True, max_length=128)
    authority = models.CharField(max_length=128, blank=True, null=True)
    description = models.TextField(blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'biodatabase'

    def get_absolute_url(self):
        return reverse('bioseq:assembly_view', args=[str(self.biodatabase_id)])