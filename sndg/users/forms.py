from django.contrib.auth import get_user_model, forms
from django.core.exceptions import ValidationError
from django.utils.translation import ugettext_lazy as _

User = get_user_model()


# class UserSignupForm(forms.Form):
#     error_message = forms.UserCreationForm.error_messages.update(
#         {"duplicate_username": _("This username has already been taken."),
#          "invalid_domain": _("Email domain not allowed"),
#          }
#     )
#
#     class Meta(forms.UserCreationForm.Meta):
#         model = User
#
#     def clean_username(self):
#         username = self.cleaned_data["username"]
#
#         data = self.cleaned_data['email']
#         domain = data.split('@')[1]
#         # if not (any([domain.endswith(x) for x in ["edu.ar", "gov.ar", "uba.ar"]])):
#         if not domain.endswith(".ar"):
#             raise ValidationError(self.error_messages["invalid_domain"])
#
#         try:
#             User.objects.get(username=username)
#         except User.DoesNotExist:
#             return username
#
#         raise ValidationError(self.error_messages["duplicate_username"])

class UserChangeForm(forms.UserChangeForm):
    class Meta(forms.UserChangeForm.Meta):
        model = User


from allauth.account.forms import SignupForm


class UserCreationForm(SignupForm):
    class Meta(forms.UserCreationForm.Meta):
        model = User

    def clean(self):
        x = super(UserCreationForm, self).clean()
        return x

    def clean_email(self):
        error_message = {"duplicate_username": _("This username has already been taken."),
                         "invalid_domain": _("Email domain not allowed"),
                         }
        error_message.update(forms.UserCreationForm.error_messages)

        username = self.cleaned_data["email"]

        data = self.cleaned_data['email']
        domain = data.split('@')[1]
        # if not (any([domain.endswith(x) for x in ["edu.ar", "gov.ar", "uba.ar"]])):
        # if not domain.endswith(".ar"):
        #     raise ValidationError(error_message["invalid_domain"])

        try:
            User.objects.get(username=username)
        except User.DoesNotExist:
            return username

        raise ValidationError(error_message["duplicate_username"])

from allauth.account.forms import LoginForm


class CustomLoginForm(LoginForm):

    def __init__(self, *args, **kwargs):
        super(CustomLoginForm, self).__init__(*args, **kwargs)
        self.fields['login'].label = _("Username")
