# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models
from django.db.models import Q

from bioseq.io.Blast import Blast
from bioresources.models.Job import Job
from bioresources.models.Barcode import Barcode

class BlastDB(models.Model):
    dbs = {"proteins": Q(biodatabase__name__endswith="_prots"),
           "rnas": Q(biodatabase__name__endswith="_rnas"),
           "barcodes": Q(biodatabase__name=Barcode.BIODBNAME),
           "structures": Q(biodatabase__name__endswith="pdb"),
           "contigs": ((~Q(biodatabase__name__endswith="_prots")) & (~Q(biodatabase__name__endswith="_rnas"))
                       & (~Q(biodatabase__name__endswith=Barcode.BIODBNAME)) & (~Q(biodatabase__name__endswith="pdb")))}

    name = models.CharField(max_length=50)
    description = models.CharField(max_length=250)
    accession = models.CharField(max_length=50)
    dbtype = models.CharField(max_length=4)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def blast(self, query_file, user, kwargs):
        Job(command=Blast.blast_cmd(), user=user)

    @staticmethod
    def init_dbs():
        BlastDB.objects.create(name="Proteins", description=_("AA sequences"), accession="proteins", dbtype="prot")
        BlastDB.objects.create(name=Barcode.BIODBNAME, description=_("Barcode Sequences"), accession="barcodes",
                               dbtype="nucl")
        BlastDB.objects.create(name="Structures", description="Protein Structure Sequences",
                               accession="structures", dbtype="prot")
        BlastDB.objects.create(name="RNAs", description="RNAs sequences", accession="rnas", dbtype="nucl")
        BlastDB.objects.create(name="Genomes", description="Genomic sequences", accession="contigs",
                               dbtype="nucl")

    class Meta:
        verbose_name_plural = _("BlastDBs")
