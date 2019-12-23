# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from sndg.users.models import User

class BatchSubmission(models.Model):

    BATCH_STATUS = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            "valid", 'invalid', "imported", "error",
        ])]
    )

    name = models.CharField(max_length=200, help_text=_('BatchName'))
    created_at = models.DateTimeField(auto_now_add=True)
    finished_at = models.DateTimeField(null=True)


    status = models.PositiveIntegerField(null=True, choices=BATCH_STATUS)


    def __str__(self):
        return "BatchSubmission(%s)" % self.name


    def __repr__(self):
        return self.__str__()


class CSVLine(models.Model):
    status = models.PositiveIntegerField(null=True, choices=BatchSubmission.BATCH_STATUS)
    line = models.TextField()
    status = models.CharField(max_length=50)
    info_status = models.CharField(max_length=250)
    batch = models.ForeignKey(BatchSubmission, related_name="csvs", on_delete=models.CASCADE)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "CSVLine('%s','%s')" % (self.db_identifier, self.process_identifier)

    def __repr__(self):
        return self.__str__()

class IdentifierLine(models.Model):
    status = models.PositiveIntegerField(null=True, choices=BatchSubmission.BATCH_STATUS)
    identifier = models.CharField(max_length=255)
    status = models.CharField(max_length=50)
    batch = models.ForeignKey(BatchSubmission, related_name="identifiers", on_delete=models.CASCADE)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "IdentifierLine('%s','%s')" % (self.db_identifier, self.process_identifier)

    def __repr__(self):
        return self.__str__()