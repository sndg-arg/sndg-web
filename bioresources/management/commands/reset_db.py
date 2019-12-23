import os
import subprocess as sp
from tqdm import tqdm
import requests
from datetime import datetime
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.core.management.base import BaseCommand

from SNDG.WebServices import download_file

from bioseq.models.Taxon import Taxon, TaxonName, TaxIdx
from bioseq.models.Ontology import Ontology
from bioseq.models.Term import Term, TermRelationship, TermDbxref, TermIdx
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

    help = 'Deletes all Resource data from the db'

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        validation = input(_("Are you  sure? If you want to continue write 'Delete'"))
        if validation == "Delete":
            self.del_db()
        else:
            self.stderr.write(_("Action canceled"))

    def del_db(self):
        from django.db import connection

        models = [Identity, ExternalId, Affiliation, ResourceRelation, ResourcePropertyValue, ResourceProperty, Person,
                  Organization,
                  Publication, Tool, Assembly, ReadsArchive, BioProject, Sample]
        pbar = tqdm(models)


        for model in pbar:
            try:
                cursor = connection.cursor()
                cursor.execute('TRUNCATE TABLE "{0}" CASCADE'.format(model._meta.db_table))
                cursor.close()
            except Exception as ex:
                pass
                print (_("Error deleting %(table)..." % {'table':model._meta.db_table}))
