# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os
import sys
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from django.conf import settings

from pdbdb.io.PDBIO import PDBIO
from bioresources.models.Job import Job
from datetime import datetime


class LoadPDBJob(Job):
    pdb = models.TextField()

    def execute(self):
        with open(self.result, "w") as stdout, open(self.dev_error, "w") as stderr:
            sys.stdout = stdout
            sys.stderr = stderr
            pdbio = PDBIO( pdbs_dir=settings.PDBSDIR ,
                 entries_path=settings.PDBSENTRIES,
                 tmp=self.job_dir())
            pdbio.init()
            pdbio.process_pdb(self.pdb)

        self.end = datetime.now()
        self.status = Job.STATUS.FINISHED
