from django.contrib import admin

from .models import Variant, Allele, Effect, Phenotype, Genotype, AntibioticResistance, Protocol, Assay, \
    ReportedAllele, GenotypeSupport, VariantCollectionSet, VariantCollectionSetAssignment, StudiedPhenotype


@admin.register(Variant)
class ToolAdmin(admin.ModelAdmin):
    search_fields = ["gene","contig"]
    # autocomplete_fields = [ "contig", "gene"]
    raw_id_fields = (
        "contig", "gene"
    )
    # list_display = ["name", "version", "url"]
    # search_fields = ["name"]


@admin.register(Allele)
class ToolAdmin(admin.ModelAdmin):
    pass


@admin.register(Effect)
class ToolAdmin(admin.ModelAdmin):
    pass


@admin.register(Phenotype)
class ToolAdmin(admin.ModelAdmin):
    pass

@admin.register(StudiedPhenotype)
class ToolAdmin(admin.ModelAdmin):
    pass

@admin.register(Genotype)
class ToolAdmin(admin.ModelAdmin):
    pass


@admin.register(AntibioticResistance)
class ToolAdmin(admin.ModelAdmin):
    list_display = ["name", "__str__"]
    raw_id_fields = (
        "antibiotic",
    )


@admin.register(Protocol)
class ToolAdmin(admin.ModelAdmin):
    pass


@admin.register(Assay)
class ToolAdmin(admin.ModelAdmin):
    raw_id_fields = (
        "result",
    )


@admin.register(ReportedAllele)
class ToolAdmin(admin.ModelAdmin):
    list_display = ["allele", "effect", "phenotype"]
    search_fields = ["phenotype__name","allele__variant__pos"]

    raw_id_fields = (
        "allele", "effect"
    )


@admin.register(GenotypeSupport)
class ToolAdmin(admin.ModelAdmin):

    raw_id_fields = (
        "reported_allele", "status","assignment",
    )


@admin.register(VariantCollectionSet)
class ToolAdmin(admin.ModelAdmin):
    pass

@admin.register(VariantCollectionSetAssignment)
class ToolAdmin(admin.ModelAdmin):
    list_display = ["variant_collection", "collection_set", "status"]

