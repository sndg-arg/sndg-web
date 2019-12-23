# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys
from django.conf import settings
from haystack import indexes

from bioseq.models.Bioentry import Bioentry
from .models.Resource import Resource
from .models.Assembly import Assembly
from .models.Organization import Organization
from .models.Person import Person
from .models.BioProject import BioProject
from .models.Publication import Publication
from .models.Structure import Structure
from .models.Tool import Tool
from .models.Expression import Expression
from .models.Barcode import Barcode
from .models.ReadsArchive import ReadsArchive

only_not_indexed = sys.argv[1] in ["rebuild_index", "update_index"]  # for solr index updating

class ResourceIndexOAI(indexes.SearchIndex, indexes.Indexable):
    """
    "item.id":1,
        "item.handle":"oai:localhost:tede/5",
        "item.lastmodified":"2018-05-05T15:09:57Z",
        "item.submitter":"submitter",
        "item.deleted":false,
        "item.public":true,
        "item.collections":["FMRP"],
        "item.communities":["com_FAMERP"],

    "metadata.dc.language":["por"],
        "metadata.dc.rights":["info:sa-repo/semantics/openAccess"],
        "metadata.dc.format":["application/pdf"],
        "metadata.dc.publisher":["Faculdade de Medicina de São José do Rio Preto",
          "Programa de Pós-Graduação em Ciências da Saúde",
          "FAMERP",
          "BR",
          "Medicina Interna; Medicina e Ciências Correlatas"],

    """
    text = indexes.CharField(document=True, use_template=True, index_fieldname=settings.HAYSTACK_DOCUMENT_FIELD)

    item_id = indexes.IntegerField(model_attr='id', index_fieldname="item.id")
    item_handle = indexes.CharField(model_attr='permalink', index_fieldname="item.handle")
    item_lastmodified = indexes.DateTimeField(model_attr='updated_at', index_fieldname="item.lastmodified")
    item_submitter = indexes.CharField(model_attr='oai_submitter', index_fieldname="item.submitter")
    item_deleted = indexes.BooleanField(model_attr='deprecated', index_fieldname="item.deleted")
    item_public = indexes.BooleanField(model_attr='oai_public', index_fieldname="item.public")
    item_collections = indexes.MultiValueField(model_attr='oai_collections', index_fieldname="item.collections")
    item_communities = indexes.MultiValueField(model_attr='oai_communities', index_fieldname="item.communities")

    item_compile = indexes.CharField(model_attr='compile', index_fieldname="item.compile")

    metadata_dc_language = indexes.MultiValueField(model_attr='metadata_dc_language',
                                                   index_fieldname="metadata.dc.language")
    metadata_dc_rights = indexes.MultiValueField(model_attr='metadata_dc_rights', index_fieldname="metadata.dc.rights")
    metadata_dc_format = indexes.MultiValueField(model_attr='metadata_dc_format', index_fieldname="metadata.dc.format")
    metadata_dc_publisher = indexes.MultiValueField(model_attr='metadata_dc_publisher',
                                                    index_fieldname="metadata.dc.publisher")

    def get_model(self):
        return Resource

    def index_queryset(self, using=None):
        """Used when the entire index for model is updated."""
        return self.get_model().objects.oai_compliant()


class PublicationIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)

    name = indexes.CharField(model_attr='name')
    description = indexes.CharField(model_attr='description')

    pub_date = indexes.DateTimeField(model_attr='date_of_publication', null=True)

    type = indexes.CharField(model_attr='type', faceted=True)

    authors = indexes.MultiValueField(model_attr='author_names', faceted=True)
    affiliations = indexes.MultiValueField(model_attr='affiliation_names', faceted=True)

    doi = indexes.CharField(model_attr='doi', null=True)
    electronic_id = indexes.CharField(model_attr='electronic_id', null=True)
    scopus_id = indexes.CharField(model_attr='scopus_id', null=True)
    issn = indexes.CharField(model_attr='issn', null=True)
    pubmed_id = indexes.CharField(model_attr='pubmed_id', null=True)

    def get_model(self):
        return Publication

    def index_queryset(self, using=None):
        return (self.get_model().objects.filter(targets__isnull=False)).distinct()


class StructureIndex(indexes.SearchIndex, indexes.Indexable):
    # TODO anotaciones GO/EC/otras...

    text = indexes.CharField(document=True, use_template=True)
    name = indexes.CharField(model_attr='name')
    description = indexes.CharField(model_attr='description')
    type = indexes.CharField(model_attr='type', faceted=True)

    pdbClass = indexes.CharField(model_attr='pdbClass', faceted=True, null=True)
    deposit_date = indexes.DateTimeField(model_attr='deposit_date', null=True)
    method = indexes.CharField(model_attr='method', faceted=True, null=True)
    org_list = indexes.CharField(model_attr='org_list', null=True)

    authors = indexes.MultiValueField(model_attr='related_author_names', faceted=True)
    affiliations = indexes.MultiValueField(model_attr='related_org_names', faceted=True)
    taxon = indexes.MultiValueField(model_attr='taxon_name', faceted=True)

    def get_model(self):
        return Structure

    def index_queryset(self, using=None):
        """Used when the entire index for model is updated."""
        return self.get_model().objects.publication_related("Argentina", only_not_indexed)


class AssemblyIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    name = indexes.CharField(model_attr='name')
    description = indexes.CharField(model_attr='description')
    type = indexes.CharField(model_attr='type', faceted=True)

    intraspecific_name = indexes.CharField(model_attr='intraspecific_name', null=True)
    species_name = indexes.CharField(model_attr='species_name', faceted=True, null=True)
    level = indexes.CharField(model_attr='level', faceted=True, null=True)
    ncbi_org = indexes.CharField(model_attr='ncbi_org', null=True)
    release_date = indexes.DateTimeField(model_attr='release_date', null=True)
    update_date = indexes.DateTimeField(model_attr='update_date', null=True)
    assembly_type = indexes.CharField(model_attr='assembly_type', faceted=True, null=True)

    authors = indexes.MultiValueField(model_attr='related_author_names', faceted=True)
    affiliations = indexes.MultiValueField(model_attr='related_org_names', faceted=True)
    taxon = indexes.MultiValueField(model_attr='taxon_name', faceted=True)

    def get_model(self):
        return Assembly

    def index_queryset(self, using=None):
        return self.get_model().objects.publication_related("Argentina", only_not_indexed)


class ExpressionIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    name = indexes.CharField(model_attr='name')
    description = indexes.CharField(model_attr='description')
    type = indexes.CharField(model_attr='type', faceted=True)

    pdat = indexes.CharField(model_attr='pdat', null=True)
    gdstype = indexes.CharField(model_attr='gdstype', faceted=True, null=True)
    submitters = indexes.CharField(model_attr='submitters', null=True)

    authors = indexes.MultiValueField(model_attr='related_author_names', faceted=True)
    affiliations = indexes.MultiValueField(model_attr='related_org_names', faceted=True)
    taxon = indexes.MultiValueField(model_attr='taxon_name', faceted=True)

    def get_model(self):
        return Expression

    def index_queryset(self, using=None):
        return self.get_model().objects.publication_related("Argentina", only_not_indexed)


class ReadsArchiveIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    name = indexes.CharField(model_attr='name')
    description = indexes.CharField(model_attr='description')
    type = indexes.CharField(model_attr='type', faceted=True)

    taxon = indexes.MultiValueField(model_attr='taxon_name', faceted=True)

    def get_model(self):
        return ReadsArchive

    def index_queryset(self, using=None):
        return self.get_model().objects.publication_related("Argentina", only_not_indexed)

class ToolIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    name = indexes.CharField(model_attr='name')
    description = indexes.CharField(model_attr='description')
    type = indexes.CharField(model_attr='type', faceted=True)

    def get_model(self):
        return Tool

    def index_queryset(self, using=None):
        return (self.get_model().objects)


class BioProjectIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    name = indexes.CharField(model_attr='name')
    description = indexes.CharField(model_attr='description')
    type = indexes.CharField(model_attr='type', faceted=True)

    sample_scope = indexes.CharField(model_attr='sample_scope', faceted=True, null=True)
    material = indexes.CharField(model_attr='material', faceted=True, null=True)
    capture = indexes.CharField(model_attr='capture', faceted=True, null=True)

    target = indexes.CharField(model_attr='target', faceted=True, null=True)
    submitters = indexes.CharField(model_attr='submitters', null=True)

    method = indexes.CharField(model_attr='method', faceted=True, null=True)

    authors = indexes.MultiValueField(model_attr='related_author_names', faceted=True)
    affiliations = indexes.MultiValueField(model_attr='related_org_names', faceted=True)
    taxon = indexes.MultiValueField(model_attr='taxon_name', faceted=True)

    def get_model(self):
        return BioProject

    def index_queryset(self, using=None):
        return self.get_model().objects.publication_related("Argentina", only_not_indexed)


class PersonIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    name = indexes.CharField(model_attr='complete_name')
    surname = indexes.CharField(model_attr='surname')

    type = indexes.CharField(model_attr='rtype', faceted=True)

    affiliations = indexes.MultiValueField(model_attr='related_org_names', faceted=True)

    def get_model(self):
        return Person

    def index_queryset(self, using=None):
        # q = { "deprecated": False }
        # if only_not_indexed:
        #     q["index_updated"] = False

        return (self.get_model().objects.filter(affiliations__resource__targets__isnull=False).distinct())


class OrganizationIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    name = indexes.CharField(model_attr='name')
    country = indexes.CharField(model_attr='country', faceted=True)
    city = indexes.CharField(model_attr='city', null=True)
    type = indexes.CharField(model_attr='rtype', faceted=True)

    def get_model(self):
        return Organization

    def index_queryset(self, using=None):
        q = {"country": "Argentina",
             "deprecated": False}
        if only_not_indexed:
            q["index_updated"] = False
        return self.get_model().objects.filter(affiliations__resource__targets__isnull=False, **q).distinct()


class BarcodeIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    name = indexes.CharField(model_attr='name')
    country = indexes.CharField(model_attr='country', faceted=True)
    subdivision = indexes.CharField(model_attr='subdivision', faceted=True, null=True)
    marker = indexes.CharField(model_attr='marker', faceted=True, null=True)
    bold_org = indexes.CharField(model_attr='bold_org', faceted=True)
    taxon = indexes.MultiValueField(model_attr='taxon_name', faceted=True)
    type = indexes.CharField(model_attr='type', faceted=True)

    def get_model(self):
        return Barcode

    def index_queryset(self, using=None):
        q = {"country": "Argentina",
             "deprecated": False}
        if only_not_indexed:
            q["index_updated"] = False
        return (self.get_model().objects.filter(**q))



# class ProteinIndex(indexes.SearchIndex, indexes.Indexable):
#     text = indexes.CharField(document=True, use_template=True)
#     locus_tag = indexes.CharField(model_attr='accession')
#     taxon = indexes.MultiValueField(model_attr='taxon__scientific_name', faceted=True)
#     description = indexes.CharField(model_attr='description')
#
#     genes = indexes.CharField(model_attr='genes', faceted=True, null=True)
#     subcellular_location = indexes.CharField(model_attr='cellular_component', faceted=True)
#     function = indexes.CharField(model_attr='molecular_function', faceted=True)
#     biological_process = indexes.CharField(model_attr='idx_biological_process', faceted=True)
#
#     # expression_condition = indexes.CharField(model_attr='ftype', faceted=True)
#     # expression_tissue = indexes.CharField(model_attr='ftype', faceted=True)
#     # interaction = indexes.CharField(model_attr='ftype', faceted=True)
#     # structure_type = indexes.CharField(model_attr='ftype', faceted=True)
#     # pathways = indexes.CharField(model_attr='ftype', faceted=True)
#
#     molecular_weight = indexes.FloatField(model_attr='molecular_weight',null=True)
#     length = indexes.IntegerField(model_attr='seq__length')
#
#     type = indexes.CharField(model_attr='ftype', faceted=True)
#
#
#     def get_model(self):
#         return Bioentry
#
#     def index_queryset(self, using=None):
#
#         return self.get_model().objects.proteins(index_updated=False)

# class SeqFeatureIndex(indexes.SearchIndex, indexes.Indexable):
#     text = indexes.CharField(document=True, use_template=True)
#     locus_tag = indexes.CharField(model_attr='locus_tag')
#     genes = indexes.CharField(model_attr='genes', faceted=True, null=True)
#     ref = indexes.CharField(model_attr='bioentry.accession')
#     description = indexes.CharField(model_attr='description')
#     taxon = indexes.MultiValueField(model_attr='bioentry.taxon.scientific_name', faceted=True)
#     type = indexes.CharField(model_attr='ftype', faceted=True)
#     feature_type = indexes.CharField(model_attr='term_type.name', faceted=True)
#     length = indexes.CharField(model_attr='ftype')
#
#     def get_model(self):
#         return Seqfeature
#
#     def index_queryset(self, using=None):
#         return self.get_model().objects.geneproducts(index_updated=False)