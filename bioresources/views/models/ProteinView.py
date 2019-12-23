# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioseq.models.Bioentry import Bioentry
from bioseq.models.Biodatabase import Biodatabase
from bioseq.models.Seqfeature import Seqfeature
from bioseq.models.Biosequence import Biosequence

from bioresources.models.Assembly import Assembly


def ProteinView(request, pk):
    be = (Bioentry.objects.select_related("biodatabase").select_related("taxon")
        .prefetch_related("dbxrefs__dbxref", "qualifiers__term", "seq",
                          "features__locations", "features__source_term",
                          "features__type_term", "features__qualifiers"))
    be = be.get(bioentry_id=pk)

    beg = Biodatabase.objects.get(name=be.biodatabase.name.replace("_prots", ""))
    assembly = Assembly.objects.get(name=beg.name)


    feature = Seqfeature.objects.seqfeature_from_locus_tag(beg.biodatabase_id, be.accession)
    feature = list(feature)[0]

    locations = list(feature.locations.all())
    start = locations[0].start_pos
    end = locations[-1].end_pos

    seq = Biosequence.objects.raw("""
    SELECT bioentry_id, version , length , alphabet ,SUBSTRING( seq,%i,%i ) seq
    FROM biosequence WHERE bioentry_id = %i ;
    """ % (start, end - start, feature.bioentry_id))[0]
    functions = {"biological_process": [], "molecular_function": [], "cellular_component": []}
    for qual in be.qualifiers.filter(term__dbxrefs__dbxref__dbname="go",
                                     term__dbxrefs__dbxref__accession="goslim_generic"):
        for dbxref in qual.term.dbxrefs.all():
            if dbxref.dbxref.accession in functions:
                functions[dbxref.dbxref.accession].append(qual.term)

    return render(request, 'resources/protein_view.html', {
        "functions": functions, "assembly":assembly,
        "object": be, "feature": feature, "taxon": assembly.ncbi_tax, "seq": seq, "start": start, "end": end,
        "sidebarleft": 1})
