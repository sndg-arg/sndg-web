import os
import subprocess as sp
from tqdm import tqdm
import requests
from datetime import datetime
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.core.management.base import BaseCommand

from bioseq.models.Taxon import Taxon, TaxonName, TaxIdx

from bioresources.models.Resource import Resource as rResource
from bioresources.models.Person import Person as rPerson
from bioresources.models.Organization import Organization as rOrganization
from bioresources.models.Affiliation import Affiliation as rAffiliation
from bioresources.models.Identity import Identity as rIdentity
from bioresources.models.ExternalId import ExternalId as rExternalId
from bioresources.models.Publication import Publication as rPublication
from bioresources.models.Expression import Expression as rExpression
from bioresources.models.Tool import Tool as rTool
from bioresources.models.Assembly import Assembly as rAssembly
from bioresources.models.ReadsArchive import ReadsArchive as rReadsArchive
from bioresources.models.BioProject import BioProject as rBioProject
from bioresources.models.Sample import Sample as rSample
from bioresources.models.Structure import Structure as rStructure
from neomodel.exceptions import AttemptedCardinalityViolation

from bioresources.models.ResourceRelation import ResourceRelation as rResourceRelation
from bioresources.models.ResourceProperty import ResourceProperty as rResourceProperty, \
    ResourcePropertyValue as rResourcePropertyValue

from django.db import transaction

# from pdbdb.models import PDB

from bioresources.graph import *


class Command(BaseCommand):
    """

    """

    help = 'Loads sndgn from the sndgr database'

    def add_arguments(self, parser):
        parser.add_argument('--remove', action='store_true')

    def handle(self, *args, **options):
        # Tax(StructuredNode):
        #     name = StringProperty()
        #     ncbi_id = StringProperty(unique_index=True)
        # for ncbi_tax_id in rResource.objects.filter(ncbi_tax__isnull=False).values_list("ncbi_tax",
        #                                                                                 flat=True).distinct():
        #     # print(ncbi_tax_id)
        #     tax = Taxon.objects.prefetch_related("names").get(ncbi_taxon_id=ncbi_tax_id)
        #     try:
        #         int(tax.scientific_name())
        #     except:
        #         print("sdflkjas")
        #         Tax.create(name=tax.scientific_name(), ncbi_id=tax.ncbi_taxon_id)

        countries = {}
        all_countries = list(rOrganization.objects.filter(
            country__isnull=False).values_list("country", flat=True).distinct())
        for country in tqdm(all_countries, file=self.stderr):
            c = Country(name=country,rid=country)
            c.save()
            countries[country] = c

        organizations = {}
        qs = rOrganization.objects.filter(country__isnull=False)
        for organization in tqdm(qs, file=self.stderr, total=qs.count()):
            org = Organization.from_resource(organization)
            org.save()
            org.location.connect(countries[organization.country])
            organizations[organization.id] = org

        # Journal(StructuredNode):
        #     name = StringProperty(unique_index=True)
        #
        #

        persons = {}
        qs = rPerson.objects.all()
        for person in tqdm(qs, file=self.stderr, total=qs.count()):
            gPerson = Person(rid=person.id,name=person.complete_name())
            gPerson.save()
            persons[person.id] = gPerson

        species = {}
        qs = rAssembly.objects.filter(
            species_name__isnull=False).values_list("species_name", flat=True).distinct()
        for species_name in tqdm(qs, file=self.stderr, total=qs.count()):
            c = Species(name=species_name)
            c.save()
            species[species_name] = c

        qs = rPublication.objects.prefetch_related("affiliations__organizations").filter(targets__isnull=False).distinct()
        for publication in tqdm(qs, file=self.stderr, total=qs.count()):
            Publication.from_resource(publication)



        qs = rExpression.objects.all()
        for expresion in tqdm(qs, file=self.stderr, total=qs.count()):
            Expression.from_resource(expresion)


        qs = rAssembly.objects.all()
        for assembly in tqdm(qs, file=self.stderr, total=qs.count()):
            Assembly.from_resource(assembly)


        qs = rStructure.objects.all()
        for structure in tqdm(qs, file=self.stderr, total=qs.count()):
            r = Structure(rid=structure.id, title=structure.name, deposit_date=structure.deposit_date,
                          method=structure.method)
            r.save()

        qs = rSample.objects.all()
        for sample in tqdm(qs, file=self.stderr, total=qs.count()):
            Sample.from_resource(sample)



        qs = rTool.objects.all()
        for tool in tqdm(qs, file=self.stderr, total=qs.count()):
            r = Tool.from_resource(tool)
            r.save()

        qs = rReadsArchive.objects.all()
        for r in tqdm(qs, file=self.stderr, total=qs.count()):
            Reads.from_resource(r)



        from neomodel import db
        resource_class = {
            rResource.RESOURCE_TYPES.STRUCTURE: Structure,
            rResource.RESOURCE_TYPES.ASSEMBLY: Assembly,
            rResource.RESOURCE_TYPES.SAMPLE: Sample,
            rResource.RESOURCE_TYPES.EXPRESSION: Expression,
            rResource.RESOURCE_TYPES.PUBLICATION: Publication,
            rResource.RESOURCE_TYPES.TOOL: Tool,
            rResource.RESOURCE_TYPES.READS: Reads,
        }

        qs = rResourceRelation.objects.all()
        for idx,rel in enumerate(tqdm(qs,total=qs.count(),file=self.stderr)):

            r = resource_class[rel.source.type]
            query = "MATCH (a:%s { rid : %i }) RETURN a" % (str(r).split(".")[-1].split("'")[0], rel.source.id)
            results, meta = db.cypher_query(query, {})
            source = r.inflate(results[0][0])

            r = resource_class[rel.target.type]
            query = "MATCH (a:%s { rid : %i }) RETURN a" % (str(r).split(".")[-1].split("'")[0], rel.target.id)
            results, meta = db.cypher_query(query, {})
            target = r.inflate(results[0][0])

            if rel.source.type == rResource.RESOURCE_TYPES.PUBLICATION:
                source.resources.connect(target)
            elif rel.source.type == rResource.RESOURCE_TYPES.SAMPLE:
                if rel.target.type == rResource.RESOURCE_TYPES.READS:
                    try:
                        target.sample.connect(source)
                    except AttemptedCardinalityViolation:
                        pass
                elif rel.target.type == rResource.RESOURCE_TYPES.ASSEMBLY:
                    try:
                        target.samples.connect(source)
                    except Exception as ex:
                        print(ex)

                else:
                    self.stderr.write(rel.role + " -> not defined")
            else:
                self.stderr.write(rel.role + " -> not defined")




        # Resource(StructuredNode):
        #     title = StringProperty(unique_index=True)
        #     description = StringProperty(unique_index=True)
        #     species = RelationshipTo('Species', 'SPECIES')
        #     tax = RelationshipTo('Tax', 'TAX')
        #

        # Barcodes(Resource):
        #     country = RelationshipTo('Country', 'COUNTRY')
        #     subdivision = StringProperty()
        #     marker = StringProperty()
        #
