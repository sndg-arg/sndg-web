import os
import subprocess as sp
from tqdm import tqdm
import requests
from datetime import datetime

from django.core.management.base import BaseCommand

from SNDG.WebServices import download_file

from bioseq.models.Taxon import Taxon, TaxonName
from bioseq.models.Ontology import Ontology
from bioseq.models.Term import Term, TermRelationship,TermDbxref
from bioseq.models.Biodatabase import Biodatabase
from bioseq.models.Bioentry import Bioentry, BioentryQualifierValue, BioentryDbxref
from bioseq.models.Biosequence import Biosequence
from bioseq.models.Seqfeature import Seqfeature, SeqfeatureDbxref, SeqfeatureQualifierValue
from bioseq.models.Location import Location
from bioseq.models.Dbxref import Dbxref, DbxrefQualifierValue

from bioresources.models.Resource import Resource
from bioresources.models.Person import Person
from bioresources.models.Organization import Organization
from bioresources.models.Affiliation import Affiliation
from bioresources.models.Identity import Identity
from bioresources.models.ExternalId import ExternalId
from bioresources.models.Publication import Publication
from bioresources.models.Tool import Tool
from bioresources.models.Assembly import Assembly
from bioresources.models.ReadsArchive import ReadsArchive
from bioresources.models.BioProject import BioProject
from bioresources.models.Sample import Sample
from bioresources.models.ResourceRelation import ResourceRelation
from bioresources.models.ResourceProperty import ResourceProperty, ResourcePropertyValue

from django.db import transaction


# from pdbdb.models import PDB


class Command(BaseCommand):
    """

    """

    help = 'Loads base taxonomy, ontologies and terms'

    def add_arguments(self, parser):
        parser.add_argument('--remove', action='store_true')

    def handle(self, *args, **options):
        with transaction.atomic():
            Ontology.load_ann_terms()
        with transaction.atomic():
            Ontology.load_go_base()

        with transaction.atomic():
            self.load_tax()



    def load_tax(self):
        pass

