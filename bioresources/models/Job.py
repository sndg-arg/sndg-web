# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from datetime import datetime
import subprocess as sp
import traceback
import os

from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from django.shortcuts import reverse
from sndg.users.models import User

from django.conf import settings
from polymorphic.models import PolymorphicModel

class Job(PolymorphicModel):
    STATUS = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            "NEW", "QUEUED", "ERROR", "RETRYING", "FINISHED", "RUNNING"
        ])]
    )

    start = models.DateTimeField(null=True)
    end = models.DateTimeField(null=True)
    result = models.TextField(null=True)
    status = models.PositiveIntegerField(default=STATUS.NEW, choices=STATUS)
    retry = models.PositiveIntegerField(default=0)
    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name="jobs", null=True)
    dev_error = models.TextField(null=True)
    result_type = models.CharField(max_length=10)

    class Meta:
        verbose_name_plural = _("Job")

    def init(self):
        assert self.id

        self.dev_error = self.job_dir() + "/err"
        self.result = self.job_dir() + "/out"
        if not os.path.exists(self.job_dir()):
            os.makedirs(self.job_dir())

    def job_dir(self):
        return settings.JOBSDIR + "/" + str(self.id) + "/"

    def queue(self):
        self.status = Job.STATUS.QUEUED
        self.start = datetime.now()

    def running(self):
        self.status = Job.STATUS.RUNNING

    def error(self):
        with open(self.dev_error, "a") as stderr:
            self.status = Job.STATUS.ERROR
            stderr.write("\n " + traceback.format_exc())

        self.end = datetime.now()

    def get_absolute_url(self):
        return reverse('bioresources:job_view', args=[str(self.id)])


class CmdJob(Job):
    command = models.TextField()

    def execute(self):
        with open(self.result, "w") as stdout, open(self.dev_error, "w") as stderr:
            sp.call(self.command, shell=True, stdout=stdout,
                    stderr=stderr)
        self.end = datetime.now()
        self.status = Job.STATUS.FINISHED
