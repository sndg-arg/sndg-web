# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os
import sys
from datetime import datetime
import gzip
import subprocess as sp
from tqdm import tqdm
from glob import glob

from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from django.conf import settings

import Bio.SeqIO as bpio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from bioseq.io.BioIO import BioIO
from bioseq.io.DB2JBrowse import DB2JBrowse
from bioresources.io.NCBISearch import NCBISearch
from bioresources.models.Job import Job
from bioresources.models.Assembly import Assembly
from bioseq.io.DB2JBrowse import DB2JBrowse


class LoadGenomeJob(Job):
    assembly = models.ForeignKey(Assembly, on_delete=models.CASCADE, related_name="jobs")

    def execute(self):
        with open(self.result, "w") as stdout, open(self.dev_error, "w") as stderr:
            self.load_genome(stderr)

            self.jbrowse_load(stderr)

        self.end = datetime.now()
        self.status = Job.STATUS.FINISHED

    def jbrowse_load(self, stderr=sys.stderr):
        io = DB2JBrowse(jbrowse_path=settings.ROOT_DIR, jbrowse_data_path=settings.SNDG_JBROWSE)
        io.stderr = stderr
        io.ovewrite = True
        io.excluded.append("gene")
        io.db2fs(self.assembly.name)


class LoadGenomeFromNCBIJob(LoadGenomeJob):

    def load_genome(self, stderr=sys.stderr):
        NCBISearch.download_assembly(self.assembly.name, workdir=self.job_dir())
        input_file = glob(self.job_dir() + "*_genomic.gbff.gz")[0]

        io = BioIO(self.assembly.name, self.assembly.ncbi_tax.ncbi_taxon_id)
        io.stderr = stderr
        grep_cmd = 'zgrep -c "FEATURES *Location/Qualifiers" "%s"' % input_file
        io.create_db()

        total = int(sp.check_output(grep_cmd, shell=True))
        with gzip.open(input_file, "rt") as h:
            io.process_record_list(bpio.parse(h, "gb"), total)


class LoadGenomeFromFileJob(LoadGenomeJob):
    filename = models.TextField()

    def load_genome(self):
        tax_id = self.assembly.ncbi_tax.ncbi_taxon_id if self.assembly.ncbi_tax else 1
        io = BioIO(self.assembly.name, tax_id)
        io.create_db()
        if self.filename.endswith("gb.gz") or self.filename.endswith("gbk.gz"):
            format = "gb"
            grep_cmd = 'zgrep -c "FEATURES *Location/Qualifiers" "%s"' % self.filename
        elif self.filename.endswith("fasta.gz") or self.filename.endswith("fna.gz"):
            format = "fasta"
            grep_cmd = 'grep -c ">" "%s"' % self.filename
        else:
            raise Exception(_("invalid file format"))
        total = int(sp.check_output(grep_cmd, shell=True))

        with gzip.open(self.filename, "rt") as h:
            io.process_record_list(bpio.parse(h, format), total)
