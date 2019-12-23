# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render
from django.db import transaction
from django.contrib.auth.decorators import login_required
from django.http import HttpResponseRedirect
from django import forms

from bioresources.models.Resource import Collaboration
from bioseq.models.Taxon import  Taxon


@login_required
def submission_start(request):
    return render(request, 'submission/submission_start.html', {})

from django_select2.forms import ModelSelect2Widget

class TaxSelect(ModelSelect2Widget):
    def label_from_instance(self,obj):
        return str(obj.scientific_name())

class TaxChoiceField(forms.ChoiceField):
    def valid_value(self, value):
        return True
    def clean(self, value):
        value = super(self.__class__,self).clean(value)
        if value:
            return Taxon.objects.get(taxon_id=value)
        return None

def form_clean_data(form):
    cleaned_data = super(form.__class__, form).clean()
    qs = form._meta.model.objects.filter(name=cleaned_data["name"])
    # if cleaned_data["ncbi_tax"]:
    #     form.data["ncbi_tax"] = Taxon.objects.get(taxon_id=cleaned_data["ncbi_tax"])
    if "pk" in form.data:
        if qs.exclude(id=form.data["pk"]).exists():
            form._errors['name'] = form._errors.get('name', [])
            form._errors['name'].append(__("%s already exists") % cleaned_data["name"])
    else:
        if qs.exists():
            form._errors['name'] = form._errors.get('name', [])
            form._errors['name'].append(__("%s already exists") % cleaned_data["name"])


def submit_model(form_class, request):
    if request.method == 'POST':
        if "pk" in request.GET:
            obj = form_class._meta.model.objects.get(id=request.GET["pk"])
            form = form_class(request.POST, instance=obj)
        else:
            form = form_class(request.POST)

        if form.is_valid():
            with transaction.atomic():
                obj = form.save()
                if not Collaboration.objects.filter(resource=obj, person=request.user.person).exists():
                    Collaboration.objects.create(resource=obj, person=request.user.person,
                                                 type=Collaboration.COLLABORATION_TYPES.owner)
            return HttpResponseRedirect(reverse("bioresources:" + obj.type_name() + "_view", args=[obj.id]))
    else:
        if "pk" in request.GET:
            resource = form_class._meta.model.objects.get(id=request.GET["pk"])
            form = form_class(instance=resource)
        else:
            form = form_class()

    data = {'form': form}
    if "pk" in request.GET:
        data["pk"] = request.GET["pk"]

    return render(request, 'submission/tool_submission.html', data)
