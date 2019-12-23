# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.contrib import admin
from django.urls import reverse
from django.utils.html import format_html

from .models.Person import Person
from .models.Affiliation import Affiliation
from .models.Organization import Organization,OrgRelationship
from .models.Resource import Resource,Collaboration
from .models.Publication import Publication
from .models.Identity import Identity
from .models.Expression import Expression
from .models.BioProject import BioProject
from .models.Structure import Structure
from .models.Sample import Sample
from .models.ReadsArchive import ReadsArchive
from .models.Assembly import Assembly
from .models.Barcode import Barcode
from .models.ResourceRelation import ResourceRelation
from .models.RKeyword import RKeyword
from .models.ExternalId import ExternalId
from .models.ResourceProperty import ResourceProperty
from .models.Job import Job
from .models.Tool import Tool



admin.site.register(Identity)
admin.site.register(RKeyword)
admin.site.register(Affiliation)
admin.site.register(ResourceProperty)


admin.site.register(ExternalId)


@admin.register(Sample)
class SampleArchiveAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]
    search_fields = ["name","description"]
    list_display = ["name", "description", "deprecated", "updated_at"]



@admin.register(Tool)
class ToolArchiveAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]
    search_fields = ["name","description"]
    list_display = ["name", "description", "deprecated", "updated_at"]

@admin.register(ReadsArchive)
class ReadsArchiveAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]
    search_fields = ["name","description"]
    list_display = ["name", "description", "deprecated", "updated_at"]

@admin.register(Resource)
class ResourceAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]
    search_fields = ["name","description"]
    list_display = ["name", "description", "deprecated", "updated_at"]

    # def link(self, obj):
    #     # https://en.proft.me/2014/10/12/reversing-admin-urls-django/
    #     return format_html('<a href="{url}?pdb_id={{pdb_id}}">Resources</a>',
    #                        pdb_id=obj.id, url=reverse('admin:bioresources_resource_changelist'))


@admin.register(ResourceRelation)
class ResourceRelationAdmin(admin.ModelAdmin):
    autocomplete_fields = ["source", "target"]


@admin.register(Expression)
class ExpressionAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]
    search_fields = ["name","description"]
    list_display = ["name", "description", "deprecated", "updated_at"]


@admin.register(BioProject)
class BioProjectAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]
    search_fields = ["name","description"]
    list_display = ["name", "description", "deprecated", "updated_at"]


@admin.register(Assembly)
class AssemblyAdmin(ResourceAdmin):
    autocomplete_fields = ["ncbi_tax"]
    search_fields = ["name","description"]
    list_display = ["name", "description", "deprecated", "updated_at"]


@admin.register(Organization)
class OrganizationAdmin(admin.ModelAdmin):
    search_fields = ["name","description","source__name"]
    autocomplete_fields = ["source"]
    list_display = ["name", "description","country","source", "deprecated", "updated_at"]

    def get_queryset(self, request):

        qs = super(OrganizationAdmin, self).get_queryset(request)

        return qs.select_related("source")


# @admin.register(OrgRelationship)
# class OrgRelationshipAdmin(admin.ModelAdmin):
#     search_fields = ["source__name","target__name","source__description","target__description"]
#     autocomplete_fields = ["source","target"]
#     list_display = ["source","target" ]
#
#     def get_queryset(self, request):
#
#         qs = super(OrgRelationshipAdmin, self).get_queryset(request)
#
#         return qs.select_related("source","target")




@admin.register(Structure)
class StructureAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]
    search_fields = ["name","description"]
    list_display = ["name", "description", "deprecated", "updated_at"]


@admin.register(Barcode)
class BarcodeAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]
    search_fields = ["name","description"]
    list_display = ["name", "description", "deprecated", "updated_at"]




@admin.register(Publication)
class PublicationAdmin(admin.ModelAdmin):
    autocomplete_fields = ["ncbi_tax"]
    list_display = ["name", "description", "links"]
    search_fields = ["name", "description"]

    def links(self, obj):
        # https://en.proft.me/2014/10/12/reversing-admin-urls-django/
        return format_html('<a href="{url}?pdb_id={{pdb_id}}">Resources</a>',
                           pdb_id=obj.id, url=reverse('admin:bioresources_resource_changelist'))

@admin.register(Job)
class JobAdmin(admin.ModelAdmin):
    autocomplete_fields = ["user"]
    list_display = ["id", "user", "status"]
    search_fields = ["name", "description"]


@admin.register(Person)
class PersonAdmin(admin.ModelAdmin):
    search_fields = [ "surname", "name"]
    list_display = [ "surname", "name"]

@admin.register(Collaboration)
class CollaborationAdmin(admin.ModelAdmin):
    autocomplete_fields = [ "person", "resource"]
    search_fields = [ "person", "resource"]
    list_display = [ "person", "resource"]