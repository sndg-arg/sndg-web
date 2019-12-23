# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.Person import Person
from bioresources.models.Resource import Resource, Collaboration
from sndg.users.models import User

from bioresources.io.GraphRepo import GraphRepo


def UserResourcesView(request,username=""):
    person = request.user.person
    collaborations = None
    graph = None
    if person:
        collaborations = set(Collaboration.objects.prefetch_related("resource").filter(person=person))
        for c in collaborations:
            c.resource.tname = Resource.RESOURCE_TYPES[c.resource.type]
        if person:
            graph, _ = GraphRepo.get_neighborhood(person.id, "Person", 1)
        else:
            graph = {}
    return render(request, 'user/user_resources.html',
                  {"user": request.user, "person": person, "collaborations": collaborations, "graph": graph})
