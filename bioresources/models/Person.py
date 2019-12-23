# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from django.urls import reverse


class Person(models.Model):

    TYPE = 20

    COLLABORATION_TYPES = Choices(
        (1, "owner", _("owner")),
        (2, "only_producer", _("only_producer")),
        (3, "only_use", _("only_use")),
        (4, "other", _("other")),
    )



    # TODO add manager to query affiliations+organizations
    surname = models.CharField(max_length=200, blank=False)
    name = models.CharField(max_length=200, default="")
    scopus_id = models.CharField(max_length=200, null=True)
    scopus_names = models.TextField(null=True, blank=True)
    email = models.EmailField(null=True, blank=True)

    deprecated = models.BooleanField(default=False)
    index_updated = models.BooleanField(default=False)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    source = models.ForeignKey("Organization", on_delete=models.SET_NULL, null=True)

    class Meta:
        verbose_name_plural = _("Persons")

    def __str__(self):
        return self.name + " " + self.surname

    def rtype(self):
        return 20  # "person"

    def get_absolute_url(self):
        return reverse('bioresources:person_view', args=[str(self.id)])

    def organizations(self):
        organizations = []
        for aff in self.affiliations.all():
            for org in aff.organizations.all():
                if org not in organizations:
                    organizations.append(org)
        return organizations

    def related_org_names(self):
        return [x.name for x in self.organizations()]

    def complete_name(self):
        return self.surname + ", " + self.name
