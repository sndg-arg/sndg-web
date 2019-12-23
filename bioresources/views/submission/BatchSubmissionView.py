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

from bioresources.models.aux.BatchSubmission import BatchSubmission


class BatchForm(forms.ModelForm):
    release_date = forms.DateField(required=True, widget=forms.SelectDateWidget(years=range(1990, datetime.now().year)))

    class Meta:
        model = BatchSubmission
        fields = ["name"]

    def __init__(self, *args, **kwargs):
        super(AssemblyForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.add_input(Submit('submit', 'Submit'))

    def clean(self):
        cleaned_data = super(ToolForm, self).clean()
        if Tool.objects.filter(name=cleaned_data["name"]).exists():
            self._errors['name'] = self._errors.get('name', [])
            self._errors['name'].append(__("%s already exists") % cleaned_data["name"])


@login_required
def BatchSubmissionView(request):
    if request.method == 'POST':
        form = BatchForm(request.POST)

        if form.is_valid():
            resource = form.save()
            return HttpResponseRedirect( reverse("Assembly_view",args=[resource.id])  )
    else:
        if "pk" in request.GET:
            resource = Assembly.objects.get(id=request.GET["pk"])
            form = AssemblyForm(instance=resource)
        else:
            form = AssemblyForm()

    return render(request, 'submission/tool_submission.html', {'form': form})

# from modeltranslation.translator import translator, TranslationOptions
# class NewsTranslationOptions(TranslationOptions):
#     fields = ('title', 'text',)
#
# translator.register(Assembly, NewsTranslationOptions)