# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models

from .Person import Person

class Identity(models.Model):
    person = models.ForeignKey(Person, on_delete=models.CASCADE)
    identifier = models.CharField(max_length=200)
    email = models.EmailField(null=True)
    url = models.URLField(null=True)
    authority = models.CharField(max_length=200)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField()
    ends = models.DateTimeField()

    class Meta:
        verbose_name_plural = _("Identities")

    def __str__(self):
        return self.identifier