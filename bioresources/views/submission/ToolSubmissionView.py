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

from bioresources.models.Tool import Tool
from bioresources.models.Resource import Collaboration
from django.db import transaction

from bioresources.views.submission import form_clean_data, submit_model

class ToolForm(forms.ModelForm):
    release_date = forms.DateField(required=True, widget=forms.SelectDateWidget(years=range(1990, datetime.now().year)))

    class Meta:
        model = Tool
        fields = ["name", "description", "url", "tool_type"]

    def __init__(self, *args, **kwargs):
        super(ToolForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.add_input(Submit('submit', 'Submit'))
        if self.instance.id:
            self.fields['name'].widget.attrs['readonly'] = True

    def clean(self):
        form_clean_data(self)


@login_required
def ToolSubmissionView(request):
    return submit_model(ToolForm, request)


