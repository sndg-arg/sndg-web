# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import redirect, reverse
from django.shortcuts import render


def sequence_view(request, pk):
    be = (Bioentry.objects.select_related("biodatabase").select_related("taxon")
        .prefetch_related("dbxrefs__dbxref", "qualifiers__term", "seq",
                          "features__locations", "features__source_term",
                          "features__type_term", "features__qualifiers"))
    be = be.get(bioentry_id=pk)

    if be.biodatabase.name.endswith("prots"):
        beg = Biodatabase.objects.get(name=be.biodatabase.name.replace("_prots", ""))
        taxon = beg.entries.first().taxon

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

    else:
        beg = be.biodatabase

    graph = entry_graph(be, beg)

    if be.biodatabase.name.endswith("prots"):
        return render(request, 'bioseq/protein_detail.html', {
            "functions": functions, "graph": graph, "accession": be.biodatabase.name.replace("_prots", ""),
            "object": be, "feature": feature, "taxon": taxon, "seq": seq, "start": start, "end": end,
            "sidebarleft": 1})
    else:
        return render(request, 'bioseq/biosequence_detail.html', {
            "object": be, "graph": graph,
            "sidebarleft": 0})  # "sidebarrigth": {"news": [{"title": "n1", "text": "lalala"}]
