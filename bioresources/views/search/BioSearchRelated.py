# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.views import View
from django.shortcuts import render

from bioresources.io.GraphRepo import GraphRepo

from bioresources.models.Resource import Resource
from bioresources.models.Person import Person
from bioresources.models.Organization import Organization
from bioresources.models.Publication import Publication
from bioresources.models.Expression import Expression
from bioresources.models.Tool import Tool
from bioresources.models.Assembly import Assembly
from bioresources.models.ReadsArchive import ReadsArchive
from bioresources.models.BioProject import BioProject
from bioresources.models.Sample import Sample
from bioresources.models.Structure import Structure


class BioSearchRelatedView(View):
    template_name = 'search/search_related.html'

    def get(self, request, rid, rtype_src, rtype_dst):
        resource_class = {
            Resource.RESOURCE_TYPES.STRUCTURE: Structure,
            Resource.RESOURCE_TYPES.ASSEMBLY: Assembly,
            Resource.RESOURCE_TYPES.SAMPLE: Sample,
            Resource.RESOURCE_TYPES.EXPRESSION: Expression,
            Resource.RESOURCE_TYPES.PUBLICATION: Publication,
            Resource.RESOURCE_TYPES.TOOL: Tool,
            Resource.RESOURCE_TYPES.READS: ReadsArchive,
            Resource.RESOURCE_TYPES.ORGANIZATION: Organization,
            Resource.RESOURCE_TYPES.PERSON: Person,
        }
        level = int(self.request.GET.get("level", "1"))

        neighborhood_ids = GraphRepo.get_neighborhood_ids(rid, rtype_src,
                                                          rtype_dst, level=1)
        if level > 1:
            neighborhood_ids += GraphRepo.get_neighborhood_ids(rid, rtype_src,
                                                              rtype_dst, level=level)

        obj = resource_class[Resource.name2code[rtype_src.upper()]].objects.get(id=rid)

        repo = resource_class[Resource.name2code[rtype_dst.upper()]]

        return render(request, self.template_name, {"rtype_dst": rtype_dst, "obj": obj,
                                                    "results": repo.objects.filter(id__in=neighborhood_ids)
                                                    })
