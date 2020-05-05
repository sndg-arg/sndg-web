# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render
from django.db import transaction

from bioresources.models.ExternalId import ExternalId
from bioresources.models.Organization import Organization
from bioresources.models.Resource import Resource, Collaboration

from bioresources.models.Person import Person
from bioresources.models.BioProject import BioProject
from bioresources.models.Sample import Sample
from bioresources.models.Expression import Expression
from bioresources.models.Publication import Publication
from bioresources.models.ResourceRelation import ResourceRelation


from bioresources.io.adapters import scopus_extended_publication


class LineResult():
    def __init__(self, line):
        self.line = line
        self.resources = []
        self.message = ""


class ResourceResult():
    def __init__(self, msg, resource):
        self.msg = msg
        self.resource = resource


resource_class = {
    Resource.RESOURCE_TYPES.PERSON: Person,
    Resource.RESOURCE_TYPES.BIOPROJECT: BioProject,
    Resource.RESOURCE_TYPES.SAMPLE: Sample,
    Resource.RESOURCE_TYPES.EXPRESSION: Expression,
    Resource.RESOURCE_TYPES.PUBLICATION: Publication,
}

from django.contrib.auth.decorators import login_required

# class RelateDoiForm(forms.Form):
#     doi = forms.CharField(max_length=100)
#
#
#     def __init__(self, *args, **kwargs):
#         super(ToolForm, self).__init__(*args, **kwargs)
#         self.helper = FormHelper()
#         self.helper.add_input(Submit('submit', 'Submit'))
#
#     def clean(self):
#         cleaned_data = super(ToolForm, self).clean()
#         qs = Tool.objects.filter(name=cleaned_data["name"])
#         # if qs.exists():
#         #     self._errors['name'] = self._errors.get('name', [])
#         #     self._errors['name'].append(__("%s already exists") % cleaned_data["name"])

from bioresources.io.scopus import scopus_client
from bioresources.io.EBISearch import EBISearch


@login_required
def relate_to_publication(request, resource_id):
    r = Resource.objects.get(id=resource_id)
    if request.method == 'POST':
        post_doi = request.POST["doi"].strip()
        qs_publication = Publication.objects.filter(doi=post_doi)
        if qs_publication.exists():
            publication = qs_publication.get()
            doi = {"title": publication.name,
                   "doi": publication.doi}
        else:
            source = "scopus"
            sds = scopus_client()
            doi = sds.doi(post_doi)
            if not doi:
                source = "ebi"
                doi = EBISearch.doi(post_doi)

        if "valid_doi" in request.POST:
            with transaction.atomic():
                if not qs_publication.exists():
                    if source == "scopus":
                        publication = scopus_extended_publication(doi["record"])
                    elif source == "ebi":
                        publication = EBISearch.save(doi)
                    else:
                        raise Exception("Invalid doi source")
                else :
                    publication = qs_publication.get()
                if not ResourceRelation.objects.filter(source=publication, target=r).count():
                    ResourceRelation.objects.create(source=publication, target=r, role="publication")

            return redirect(r.get_absolute_url())
        else:
            data = {'resource': r}
            if doi:
                data['valid_doi'] = "yes"
                data['doi'] = post_doi
                data['msg'] = doi["title"]
            else:
                data['msg'] = "invalid doi:" + post_doi

            return render(request, 'submission/relate_to_publication.html', data)

    return render(request, 'submission/relate_to_publication.html', {'resource': r})


@login_required
def mark_to_relate(request, resource_id):
    request.session["relate_with"] = resource_id
    return redirect(request.POST["next"])


@login_required
def claim_identity(request, person_id):
    p = Person.objects.get(id=person_id)
    if request.method == 'POST':
        request.user.person = p
        request.user.save()
        return redirect(reverse('bioresources:user_resources'))
    else:
        return render(request, 'submission/claim_identity.html', {'person': p})


@login_required
def claim_resource(request, resource_id):
    from bioresources.graph import connect_nodes
    r = Resource.objects.get(id=resource_id)

    if request.method == 'POST':
        person = request.user.person
        c = Collaboration(person=person, resource=r,
                          type=request.POST["relation"])
        c.save()
        person.type = Resource.RESOURCE_TYPES.PERSON
        x = connect_nodes(person, r, reltype=str(Collaboration.COLLABORATION_TYPES[int(c.type)]))
        return redirect(r.get_absolute_url())
    else:
        return render(request, 'submission/claim_resource.html', {
            "collaboration_types": {x[0]: str(x[1]) for x in
                                    Collaboration.COLLABORATION_TYPES},
            'resource': r})


@login_required
def SubmissionRelatedView(request, src_id, dst_id):
    from bioresources.graph import connect_nodes
    r1 = Resource.objects.get(id=src_id)
    r2 = Resource.objects.get(id=dst_id)

    if request.method == 'POST':
        rr = ResourceRelation(source=r1, target=r2, role="uses")
        rr.save()
        connect_nodes(r1, r2)

        return redirect(r2.get_absolute_url())
    else:

        return render(request, 'submission/submission_related.html', {
            'resource1': r1, "resource2": r2})
