from django.shortcuts import redirect,reverse
from django.http import HttpResponseNotFound
from bioseq.models.Bioentry import Bioentry


def search_redirect(request,acctype,acc):
    if acctype == "locus_tag":
        be = Bioentry.objects.get(accession=acc)
        return redirect(reverse("bioresources:protein_view",args=[be.bioentry_id]))
    else:
        return HttpResponseNotFound("%s -> %s" % (acctype,acc))