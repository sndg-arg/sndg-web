import datetime

from captcha.fields import CaptchaField
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit
from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from haystack.forms import SearchForm

from .models import Resource, Assembly






class PublicationForm(forms.Form):
    name = forms.CharField(max_length=100)
    type = forms.CharField(max_length=100)
    start_date = forms.DateField()
    end_date = forms.DateField()

    def clean(self):
        cleaned_data = super(PublicationForm, self).clean()
        start_date = cleaned_data.get("start_date")
        end_date = cleaned_data.get("end_date")

        if start_date and end_date and (start_date > end_date):
            self._errors['start_date'] = self._errors.get('start_date', [])
            self._errors['start_date'].append("Start date must be before end date.")

        return cleaned_data


class SignUpForm(UserCreationForm):
    first_name = forms.CharField(max_length=30, required=False, help_text='Optional.')
    last_name = forms.CharField(max_length=30, required=False, help_text='Optional.')
    email = forms.EmailField(max_length=254, help_text='Required. Inform a valid email address.')

    class Meta:
        model = User
        fields = ('username', 'first_name', 'last_name', 'email', 'password1', 'password2',)


class AllauthSignupForm(forms.Form):
    captcha = CaptchaField()

    def signup(self, request, user):
        """ Required, or else it throws deprecation warnings """
        pass


