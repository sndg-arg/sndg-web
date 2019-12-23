from Bio.Seq import Seq
from biosql.models import Biodatabase, Bioentry, Seqfeature, Term, Ontology
from django.db import models
from django.db.models import Q
from django.utils.translation import gettext as _
from model_utils import Choices

from .managers import VariantannotationManager, ReportedAlleleManager, VariantassignmentManager


class Allele(models.Model):
    id = models.AutoField(primary_key=True)
    variant_fk = models.ForeignKey('Variant', models.DO_NOTHING, db_column='variant_fk', related_name="alleles")
    alt = models.CharField(max_length=255)
    main_effect_fk = models.ForeignKey('Effect', models.DO_NOTHING, db_column='main_effect_fk')
    hgvs_c = models.CharField(max_length=255)

    def check_effect(self, seq_str):
        mut = Seq(seq_str[:self.variant_fk.pos] + self.alt + seq_str[self.variant_fk.pos + len(self.alt) - len(
            self.variant_fk.ref):])
        loc = self.main_effect_fk.gene.locations.all()[0]
        start, end = loc.start_pos, loc.end_pos
        mut = mut[start:end]
        if loc.strand == -1:
            mut = mut.reverse_complement()

        return self.main_effect_fk.aa_alt == str(mut.translate())[self.main_effect_fk.aa_pos + 1]

    class Meta:
        unique_together = (('variant_fk', 'alt'),)

    def __str__(self):
        return str(self.variant_fk) + "->" + self.alt


class AlleleEffect(models.Model):
    allele_fk = models.ForeignKey(Allele, models.CASCADE, db_column='allele_fk', related_name="effects")
    effect_fk = models.ForeignKey('Effect', models.CASCADE, db_column='effect_fk', related_name="alleles")


class Effect(models.Model):
    id = models.AutoField(primary_key=True)

    transcript = models.CharField(max_length=255, blank=True, db_index=True)
    gene = models.ForeignKey(Seqfeature, models.CASCADE, related_name="effects")
    variant_type = models.CharField(max_length=255)

    predicted_impact = models.CharField(max_length=255, blank=True, null=True)
    aa_ref = models.CharField(max_length=255, blank=True, null=True)
    aa_pos = models.IntegerField(blank=True, null=True)
    aa_alt = models.CharField(max_length=255, blank=True, null=True)
    hgvs_p = models.CharField(max_length=255, blank=True, null=True)

    ref_organism = models.ForeignKey(Biodatabase, models.CASCADE, "variant_effects", db_index=True)

    def __str__(self):
        return self.transcript + " " + self.variant_type + " " + self.hgvs_p


class Variant(models.Model):
    id = models.AutoField(primary_key=True)
    pos = models.IntegerField()
    # gene = models.CharField(max_length=255, blank=True, null=True)
    gene = models.ForeignKey(Seqfeature, models.CASCADE, related_name="variants")
    gene_pos = models.IntegerField(blank=True)

    description = models.CharField(max_length=255, blank=True, null=True)

    contig = models.ForeignKey(Bioentry, models.CASCADE, related_name="variants")
    ref = models.CharField(max_length=255)

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return self.contig.name + " " + self.ref + str(self.pos)

    class Meta:
        unique_together = (('contig', 'pos', 'ref'),)
        index_together = [
            ["contig", "gene"],
        ]


class Variantannotation(models.Model):
    DP = "DP"

    id = models.AutoField(primary_key=True)
    assignment_fk = models.ForeignKey('Variantassignment', models.DO_NOTHING, db_column='assignment_fk',
                                      related_name="annotations")
    source_type = models.CharField(max_length=255)
    source = models.CharField(max_length=255)
    prop = models.CharField(max_length=255)
    value = models.CharField(max_length=255)
    description = models.CharField(max_length=255, blank=True, null=True)

    objects = VariantannotationManager()

    def __str__(self):
        return self.prop + "=" + self.value


class Variantassignment(models.Model):
    id = models.AutoField(primary_key=True)
    variant_collection_fk = models.ForeignKey('Variantcollection', models.CASCADE, related_name="assignments",
                                              db_column='variant_collection_fk')
    variant_fk = models.ForeignKey(Variant, models.CASCADE, db_column='variant_fk', related_name="assignments")
    allele_fk = models.ForeignKey(Allele, models.CASCADE, db_column='allele_fk', related_name="assignments")

    objects = VariantassignmentManager()

    class Meta:
        unique_together = (('variant_collection_fk', 'variant_fk', 'allele_fk'),)


class Variantcollection(models.Model):
    id = models.AutoField(primary_key=True)
    ref_organism = models.ForeignKey(Biodatabase, models.CASCADE, "strains")
    sample = models.CharField(max_length=255)
    description = models.CharField(max_length=255, default="")
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return self.sample

    def __repr__(self):
        return str(self)

    class Meta:
        unique_together = (('ref_organism', 'sample'),)


class Phenotype(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=255)
    description = models.TextField()

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)


class AntibioticResistance(Phenotype):
    antibiotic = models.ForeignKey(Term, models.DO_NOTHING)

    def __str__(self):
        return self.antibiotic.identifier

    def __repr__(self):
        return str(self)

    def process_variant_collection(self, genotype, variant_collection, depth=30):
        """
        processes a variant collection and creates a GenotypeSupport if needed
        :param vc: variant collection
        :return:
        """

        support = []
        status_ontology = Ontology.objects.get(name=GenotypeSupport.STATUS_ONTOLOGY)
        status_conclusive = Term.objects.get(ontology=status_ontology, name="Conclusive")
        status_possible = Term.objects.get(ontology=status_ontology, name="Possible")
        status_hint = Term.objects.get(ontology=status_ontology, name="Hint")

        # low_conf_alleles = [x.assignment_fk.allele_fk for x in
        #                     Variantannotation.objects.dp_filter(variant_collection, depth).prefetch_related(
        #                         "assignment_fk__allele_fk")]
        low_conf_alleles = []
        exact_reported = list(ReportedAllele.objects.exact_reported(self, variant_collection).exclude(
            allele__in=low_conf_alleles))

        support += self.get_genotype_support(exact_reported, variant_collection, status_conclusive, genotype)

        pos_reported = [x for x in ReportedAllele.objects.pos_reported(self, variant_collection).exclude(
            allele__in=low_conf_alleles) if x not in exact_reported]

        support += self.get_genotype_support(pos_reported, variant_collection, status_possible, genotype)

        reported_genes = self.reported_genes()
        gene_variants = Variantassignment.objects.reported_positions(variant_collection,
                                                                     reported_genes).exclude(
            allele_fk__in=low_conf_alleles)

        for assignment in gene_variants:
            if (assignment.allele_fk not in [r.allele for r in exact_reported]) and (
                    assignment.allele_fk not in [r.allele for r in pos_reported]):

                if (assignment.allele_fk.main_effect_fk and
                        assignment.allele_fk.main_effect_fk.variant_type == "frameshift_variant"):
                    gene = [x for x in reported_genes if x == assignment.allele_fk.variant_fk.gene][0]
                    variant_fs = [x for x in gene.effects.all() if "frameshift_variant" in x.variant_type]

                    if variant_fs:
                        gs = GenotypeSupport(genotype=genotype, status=status_conclusive, assignment=assignment)
                    else:
                        gs = GenotypeSupport(genotype=genotype, status=status_possible, assignment=assignment)

                else:
                    gs = GenotypeSupport(genotype=genotype, status=status_hint, assignment=assignment)

                support.append(gs)
        return support

    def get_genotype_support(self, reported_qs, variant_collection, status_term, genotype):
        support = []
        for reported in reported_qs:
            if reported.allele:

                assignment = Variantassignment.objects.get(variant_collection_fk=variant_collection,
                                                           variant_fk=reported.allele.variant_fk)
            else:
                alleles = [x.allele_fk for x in reported.effect.alleles.all()]
                assignment = Variantassignment.objects.get(variant_collection_fk=variant_collection,
                                                           allele_fk__in=alleles)

            gs = GenotypeSupport(genotype=genotype, status=status_term, reported_allele=reported, assignment=assignment)
            support.append(gs)
        return support

    def reported_genes(self):

        return Seqfeature.objects.prefetch_related("qualifiers__term").filter(
            Q(variants__alleles__main_effect_fk__reported__phenotype=self) |
            Q(variants__alleles__reported__phenotype=self)).distinct()


class Protocol(models.Model):
    id = models.AutoField(primary_key=True)
    phenotype = models.ForeignKey(Phenotype, models.CASCADE, related_name="protocols")
    name = models.CharField(max_length=255)

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)


class Assay(models.Model):
    ASSAY_ONTOLOGY = "Assay Results"

    id = models.AutoField(primary_key=True)
    protocol = models.ForeignKey(Protocol, models.CASCADE, related_name="assays")
    variant_collection = models.ForeignKey(Variantcollection, models.CASCADE, related_name="assays")
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    result = models.ForeignKey(Term, models.DO_NOTHING)

    def __str__(self):
        return "%s %s %s" % (self.variant_collection.sample, self.protocol.phenotype.name, self.result.identifier)

    def __repr__(self):
        return str(self)


class Genotype(models.Model):
    id = models.AutoField(primary_key=True)
    phenotype = models.ForeignKey(Phenotype, models.CASCADE, related_name="genotypes")
    variant_collection = models.ForeignKey(Variantcollection, models.CASCADE,related_name="genotypes")

    def __str__(self):
        return "%s %s" % (self.phenotype.name, self.variant_collection.sample)


class ReportedAllele(models.Model):
    id = models.AutoField(primary_key=True)
    allele = models.ForeignKey(Allele, models.CASCADE, null=True, related_name="reported")
    effect = models.ForeignKey(Effect, models.CASCADE, related_name="reported")
    phenotype = models.ForeignKey(Phenotype, models.CASCADE, "reported")
    reported_in = models.CharField(max_length=255)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    objects = ReportedAlleleManager()

    class Meta:
        unique_together = (('effect', 'phenotype'),)

    def __str__(self):
        return (("%s %s %s" % (self.allele.hgvs_c, self.phenotype.name, self.effect.variant_type))
        if self.allele else ("%s %s" % (str(self.effect), self.phenotype.name)))


class GenotypeSupport(models.Model):
    STATUS_ONTOLOGY = "Genotype Support Status"

    id = models.AutoField(primary_key=True)
    genotype = models.ForeignKey(Genotype, models.CASCADE, related_name="supportedby")
    status = models.ForeignKey(Term, models.DO_NOTHING)
    reported_allele = models.ForeignKey(ReportedAllele, models.DO_NOTHING, null=True)
    assignment = models.ForeignKey(Variantassignment, models.DO_NOTHING)

    def __str__(self):
        return "%s ->  %s (%s)" % (self.genotype.variant_collection.sample,
                                   str(self.reported_allele), self.status)


class VariantCollectionSet(models.Model):
    STAT_METRICS = ["sensitivity", "specificity", "TP", "TN", "FN", "FP"]

    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=255)
    description = models.TextField(default="")

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)

    def samples(self):
        return [x.variant_collection for x in self.assignments.all() if
                x.status == VariantCollectionSetAssignment.VCSTATUS.active]

    def phenotypes(self):
        return [x.phenotype for x in self.studied_phenotypes.all() if
                x.status == StudiedPhenotype.PSTATUS.active]





class VariantCollectionSetAssignment(models.Model):
    VCSTATUS = Choices((1, "active", _("+")),
                       (2, "inactive", _("-")),
                       )
    variant_collection = models.ForeignKey(Variantcollection, models.CASCADE, related_name="collection_sets")
    collection_set = models.ForeignKey(VariantCollectionSet, models.CASCADE, related_name="assignments")
    status = models.PositiveIntegerField(choices=VCSTATUS, default=VCSTATUS.active)


class StudiedPhenotype(models.Model):
    PSTATUS = Choices((1, "active", _("+")),
                      (2, "inactive", _("-")),
                      )
    phenotype = models.ForeignKey(Phenotype, models.CASCADE, related_name="studiedin")
    collection_set = models.ForeignKey(VariantCollectionSet, models.CASCADE, related_name="studied_phenotypes")
    status = models.PositiveIntegerField(choices=PSTATUS, default=PSTATUS.active)

    def __str__(self):
        return "%s %s %s" % (self.phenotype.name, self.collection_set.name, StudiedPhenotype.PSTATUS[self.status])

    def __repr__(self):
        return str(self)
