import json
import math
from collections import Counter

from biosql.models import Term
from django.db.models import Q, Count
from django.http import HttpResponseBadRequest
from django.shortcuts import render

from .io.PhenoGenoTable import PhenoGenoTable
from .models import Variantcollection, Phenotype, Genotype, Assay, Allele, Effect, Variant, \
    Variantassignment, VariantCollectionSet, GenotypeSupport


class Page:
    def __init__(self):
        self.count = None
        self.size = None
        self.offset = None
        self.page = None
        self.number = None
        self.previous_page_number = None
        self.next_page_number = None
        self.has_previous = None
        self.has_next = None
        self.num_pages = None
        self.paginator = {}


def variant_collection(request, pk):
    vc = Variantcollection.objects.get(id=pk)

    # Phenotypes / Genotypes
    pheno_table = []
    reported_alleles = {}
    for p in Phenotype.objects.filter(Q(genotypes__variant_collection=vc) |
                                      Q(protocols__assays__variant_collection=vc)):
        pheno = {"name": p.name}
        gquery = Genotype.objects.prefetch_related("supportedby__reported_allele__allele").filter(variant_collection=vc,
                                                                                                  phenotype=p)
        if gquery.exists():
            g = gquery.get()
            pheno["genotype"] = list(g.supportedby.all())
            reported_alleles.update(
                {x.reported_allele.allele: x.reported_allele for x in pheno["genotype"] if x.reported_allele})

        else:
            pheno["genotype"] = []
        aquery = Assay.objects.filter(variant_collection=vc, protocol__phenotype=p)
        if aquery.exists():
            a = aquery.get()
            pheno["phenotype"] = a.result.identifier
        else:
            pheno["phenotype"] = "?"

        if aquery.exists() or gquery.exists():
            pheno_table.append(pheno)

    # summaries
    Effect.objects.annotate(tcount=Count("variant_type")).filter(allele__assignments__variant_collection_fk=vc)

    # List
    alleles = Allele.objects.annotate().prefetch_related("variant_fk__gene__qualifiers", "main_effect_fk").filter(
        assignments__variant_collection_fk=vc)
    alleles = alleles[0:20]
    for a in alleles:
        if a in reported_alleles:
            a.reported_a = reported_alleles[a]
        else:
            a.reported_a = ""

    return render(request, 'variant_collection.html', {
        "variant_collection": vc, "alleles": alleles,
        "pheno_table": pheno_table})


def variant(request, pk):
    variant = Variant.objects.get(id=pk)

    if "variant_collection" in request.GET:
        vc = Variantcollection.objects.get(id=int(request.GET["variant_collection"]))
        va = list(Variantassignment.objects.filter(variant_collection_fk=vc, variant_fk=variant))
        return render(request, 'variant.html', {
            "assignments": va, "hgvs": "%s:%i-%s" % (variant.contig.name, variant.pos, variant.ref)})
    elif "collection_set" in request.GET:
        return render(request, 'variant.html', {
            "allele": allele})
    else:
        raise HttpResponseBadRequest("either variant_collection or collection_set are required")


def collection_set(request, pk):
    vcs = VariantCollectionSet.objects.prefetch_related(
        "assignments__variant_collection__genotypes__supportedby__reported_allele__allele",
        "assignments__variant_collection__assays__protocol__phenotype", "studied_phenotypes"
    ).get(id=pk)

    pgt = PhenoGenoTable()
    sample_pheno = pgt.sample_pheno(vcs)

    pos = Term.objects.get(ontology__name=Assay.ASSAY_ONTOLOGY, identifier="Positive")
    neg = Term.objects.get(ontology__name=Assay.ASSAY_ONTOLOGY, identifier="Negative")

    conclusive = Term.objects.get(ontology__name=GenotypeSupport.STATUS_ONTOLOGY, identifier="Conclusive")
    possible = Term.objects.get(ontology__name=GenotypeSupport.STATUS_ONTOLOGY, identifier="Possible")
    hint = Term.objects.get(ontology__name=GenotypeSupport.STATUS_ONTOLOGY, identifier="Hint")

    stats, stats_genetic, pheno_sample = pgt.pheno_sample_stats(vcs, pos, neg, conclusive, possible, hint)
    _, stats_genetic_conclusive, _ = pgt.pheno_sample_stats(vcs, pos, neg, conclusive, possible, hint,
                                                            lambda g: g.status == conclusive)
    for p, gts in list(stats_genetic_conclusive.items()):
        for gt, ms in list(gts.items()):
            if len(ms["TP"]) == 0 and len(ms["FP"]) == 0:
                del stats_genetic_conclusive[p][gt]
        stats_genetic_conclusive[p] = dict(stats_genetic_conclusive[p])
    stats_genetic_conclusive = dict(stats_genetic_conclusive)

    # # summaries
    # Effect.objects.annotate(tcount=Count("variant_type")).filter(allele__assignments__variant_collection_fk=vc)

    # List

    queryset = Variant.objects.prefetch_related(
        "assignments__allele_fk__reported", "contig", "assignments__annotations",
        "gene__qualifiers", "assignments__allele_fk__main_effect_fk").filter(
        # assignments__annotations__prop="DP",assignments__annotations__value__gt=50,
        assignments__variant_collection_fk__in=vcs.samples()).annotate(
        alleles_count=Count('assignments__allele_fk__alt'))

    # queryset = Variantassignment.objects.prefetch_related("allele_fk__reported","variant_fk__contig",
    #     "annotations","variant_fk__gene__qualifiers",
    #     "allele_fk__main_effect_fk").filter(variant_collection_fk__in=vcs.samples())

    if "only_diff" in request.GET:
        queryset = queryset.filter(alleles_count__lt=len(vcs.samples()))
    queryset = queryset.distinct()
    page = page_from_request(queryset.count(), request)
    queryset = queryset[page.offset:page.offset + page.size]

    variant_dict = []
    for v in queryset:
        alleles_dict = {"pos": v.pos, "chrom": v.contig.name, "ref": v.ref, "id": v.id,
                        "gene": v.gene.qualifiers_dict()["locus_tag"], "gene_pos": v.gene_pos, "samples": []}

        for s in vcs.samples():
            assignment = [va for va in v.assignments.all() if va.variant_collection_fk == s]
            if assignment:
                va = assignment[0].allele_fk
                sample_dict = {"alt": va.alt, "effect": va.main_effect_fk,
                               "annotations": {x.prop:x.value for x in assignment[0].annotations.all()}}
            else:
                sample_dict = {"alt": v.ref, "effect": None, "annotations": []}
            sample_dict["name"] = s.sample
            alleles_dict["samples"].append(sample_dict)

        alleles_dict.update(Counter([x["alt"] for x in alleles_dict["samples"]]))

        variant_dict.append(alleles_dict)

    p_total = stats["total"]
    del stats["total"]

    query = ""
    if "format" in request.GET and request.GET["format"] == "json":
        data = {
            "variants": variant_dict
        }
        from django.http import HttpResponse
        return HttpResponse(json.dumps(data, indent=4, sort_keys=True,default=lambda obj:obj.__dict__))
    else:
        return render(request, 'variant_collection_set.html', {
            "variant_collection_set": vcs, "variants": variant_dict,
            "pheno_table": sample_pheno, "samples": vcs.samples(), "stats_genetic": stats_genetic,
            "stats": stats, "pheno_sample": pheno_sample, "p_total": p_total, "query": query, "page_obj": page,
            "phenotypes": vcs.phenotypes(), "metrics": VariantCollectionSet.STAT_METRICS,
            "stats_genetic_conclusive": stats_genetic_conclusive
        })


def page_from_request(count, request):
    page_size = int(request.GET.get("page_size", "20"))
    if page_size > 500:
        page_size = 500
    page_num = int(request.GET.get("page", "1"))
    page_offset = (page_num - 1) * page_size

    prange = range(1, math.ceil(count / page_size))
    number = math.floor(page_offset / page_size) + 1
    page = Page()
    page.count = count
    page.size = page_size
    page.offset = page_offset
    page.page = page_num
    page.number = number
    page.previous_page_number = number - 1
    page.next_page_number = number + 1
    page.has_previous = page_num > 1
    page.has_next = page_num < math.ceil(count / page_size)
    page.num_pages = math.ceil(count / page_size)
    page.paginator = {"show_pages": (prange[max(
        number - 3, 0):number + 2]),
                      "num_pages": math.ceil(count / page_size), "count": count}

    return page
