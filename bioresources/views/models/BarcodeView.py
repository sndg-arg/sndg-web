# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioresources.models.Barcode import Barcode
from bioresources.io.GraphRepo import GraphRepo

def barcode(request, pk):
    barcode = Barcode.objects.get(id=pk)
    graph, related_resources = GraphRepo.get_neighborhood(pk, "Barcode",level=1)
    external_url ="http://www.boldsystems.org/index.php/Public_RecordView?processid=" + barcode.name
    return render(request, 'resources/barcode.html', { "external_url":external_url,
        "graph": graph, "related_resources": related_resources, "pk": pk, "rtype_src": "Barcode",
        "barcode": barcode,     "sidebarleft": 1, })
