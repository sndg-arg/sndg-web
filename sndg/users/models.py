from django.contrib.auth.models import AbstractUser
from django.db import models
from django.urls import reverse
from django.utils.translation import ugettext_lazy as _
from bioresources.models.Person import Person
from bioresources.models.Resource import Collaboration


class User(AbstractUser):
    person = models.OneToOneField(Person, on_delete=models.SET_NULL, related_name="user", null=True)
    conicet_id = models.IntegerField(null=True)

    def get_collaboration(self, resource):
        if bool(self.person):
            collaboration = list(Collaboration.objects.filter(person=self.person, resource=resource))
            if collaboration:
                return collaboration[0]
        else:
            return None

    def validated(self):
        return self.conicet_id != None

    def get_absolute_url(self):
        return reverse("users:detail", kwargs={"username": self.username})


class UserPerson(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name="persons")
    person = models.ForeignKey(Person, on_delete=models.CASCADE, related_name="users")
