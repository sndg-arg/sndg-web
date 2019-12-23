# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.views.generic import TemplateView

from bioseq.models.Taxon import Taxon

class TaxView(TemplateView):
    model = Taxon
    template_name = "resources/taxon_view.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['object'] = Taxon.objects.prefetch_related("names").get(ncbi_taxon_id=self.kwargs["pk"])
        return context