# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models

from .Publication import Publication
from .Resource import Resource
from .Person import Person
from .Organization import Organization

class Affiliation(models.Model):
    resource = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="affiliations")
    author = models.ForeignKey(Person, on_delete=models.CASCADE, related_name="affiliations")
    organizations = models.ManyToManyField(Organization, related_name="affiliations")

    class Meta:
        verbose_name_plural = _("Affiliations")

    def __str__(self):
        # resource = models.ForeignKey(Resource, on_delete=models.CASCADE, related_name="affiliations")
        orgs = "-".join( [x.name for x in self.organizations.all()])
        author = str(self.author)
        r = self.resource.name
        return ("Affiliation: "  + author + "|" + r + "|" + orgs)