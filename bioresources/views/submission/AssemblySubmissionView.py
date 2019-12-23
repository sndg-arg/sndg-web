# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from datetime import datetime

from django.shortcuts import redirect, reverse
from django.shortcuts import render
from django import forms

from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit
from django.http import HttpResponseRedirect
from django.contrib.auth.decorators import login_required

from bioresources.models.Assembly import Assembly
from bioresources.views.submission import form_clean_data, submit_model

from . import TaxChoiceField, TaxSelect

from bioseq.models.Taxon import Taxon


class TaxChoiceField2(forms.ChoiceField):
    def valid_value(self, value):
        return True
    def clean(self, value):
        value = super(self.__class__,self).clean(value)
        if value:
            return Taxon.objects.get(taxon_id=value).scientific_name()
        return None

class AssemblyForm(forms.ModelForm):
    release_date = forms.DateField(required=True, widget=forms.SelectDateWidget(years=range(1990, datetime.now().year)))

    ncbi_tax = TaxChoiceField(
        widget=TaxSelect(
            model=Taxon,
            search_fields=['names__name__icontains']
        ),required=False
    )

    species_name = TaxChoiceField2(
        widget=TaxSelect(
            queryset=Taxon.objects.filter(node_rank="species"),
            model=Taxon,
            search_fields=['names__name__icontains']
        ),required=False
    )


    class Meta:
        model = Assembly
        fields = ["name", "description", "intraspecific_name", "species_name", "level", "ncbi_tax"]

    def __init__(self, *args, **kwargs):
        super(AssemblyForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.add_input(Submit('submit', 'Submit'))
        if self.instance.id:
            self.fields['name'].widget.attrs['readonly'] = True

    def clean(self):
        form_clean_data(self)


@login_required
def AssemblySubmissionView(request):
    return submit_model(AssemblyForm, request)

# from modeltranslation.translator import translator, TranslationOptions
# class NewsTranslationOptions(TranslationOptions):
#     fields = ('title', 'text',)
#
# translator.register(Assembly, NewsTranslationOptions)
