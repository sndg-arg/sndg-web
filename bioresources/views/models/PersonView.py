# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.Person import Person
from bioresources.io.GraphRepo import GraphRepo





def person(request, pk):
    person = Person.objects.get(id=pk)
    graph, related_resources = GraphRepo.get_neighborhood(pk, "Person")
    return render(request, 'resources/person.html', {"pk": person.id, "related_resources": related_resources,
                                                     "person": person, "graph": graph,
                                                     "sidebarleft": 1, })
