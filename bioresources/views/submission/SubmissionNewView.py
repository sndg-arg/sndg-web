# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from datetime import datetime

from django.shortcuts import redirect, reverse
from django.shortcuts import render
from django import forms

from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit,Field


from bioresources.models.Resource import UnprocessedResource,Collaboration
from bioseq.models.Taxon import Taxon
from django.contrib.auth.decorators import login_required

class UnprocessedResourceForm(forms.ModelForm):
    #     # name = forms.CharField(max_length=350,required=True)
    #     # description = forms.CharField(widget=forms.Textarea,required=False)
    #     # intraspecific_name = forms.CharField(max_length=250, required=False)
    #     # species_name = forms.CharField(max_length=200, required=False)
    #     # level = forms.CharField(max_length=50, required=True)
    #     # # ncbi_org = forms.CharField(max_length=200, null=True)
    # release_date = forms.DateField(required=True, widget=forms.SelectDateWidget(years=range(1990, datetime.now().year)))

    #     # # update_date = forms.DateField(null=True)
    #     # assembly_type = forms.ChoiceField( required=True, choices=(
    #     #     ("haploid", "haploid"), ("diploid", "diploid"), ("other", "other")))
    relation = forms.ChoiceField( required=True, choices=Collaboration.COLLABORATION_TYPES )
    ncbi_tax = forms.IntegerField(required=False)
    class Meta:
        model = UnprocessedResource
        fields = ["name", "description", "future_type","ncbi_tax"]  # ,"ncbi_tax"

    def __init__(self, *args, **kwargs):
        super(UnprocessedResourceForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.add_input(Submit('submit', 'Submit'))

        self.fields["future_type"].choices = [x for x in list(self.fields["future_type"].choices) if
                                               (not x[0]) or (x[0] < 11)]

    def clean(self):
        cleaned_data = super(UnprocessedResourceForm, self).clean()
        qs = Taxon.objects.filter(ncbi_taxon_id=cleaned_data["ncbi_tax"])
        if qs.count():
            cleaned_data["ncbi_tax"] = qs.get()
        else:
            self._errors['ncbi_tax'] = self._errors.get('ncbi_tax', [])
            self._errors['ncbi_tax'].append(_("%s does not exists") % cleaned_data["ncbi_tax"])


        if UnprocessedResource.objects.filter(name=cleaned_data["name"]).exists():
            self._errors['name'] = self._errors.get('name', [])
            self._errors['name'].append(_("%s already exists") % cleaned_data["name"])



@login_required
def SubmissionNewView(request):
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = UnprocessedResourceForm(request.POST)

        if form.is_valid():
            record = form.save()
            relation = form.cleaned_data["relation"]
            Collaboration(resource=record, person=request.user.person, type=relation).save(force_insert=True)
            return HttpResponseRedirect(reverse("bioresources:user_resources"))


    else:
        form = UnprocessedResourceForm()

    return render(request, 'submission/submission_new.html', {'form': form})



