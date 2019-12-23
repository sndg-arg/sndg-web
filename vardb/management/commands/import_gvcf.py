import gzip
import logging
import subprocess
from functools import reduce

import hgvs.parser
import vcf
from biosql.models import Bioentry, Biodatabase, Seqfeature
from django.db import transaction
from django.db.models import Q
from django_tqdm import BaseCommand
from tqdm import tqdm

from vardb.models import Variant, Allele, Variantannotation, Variantassignment, Variantcollection, Effect, AlleleEffect

_log = logging.getLogger(__name__)


# log = logging.getLogger('django.db.backends')
# log.setLevel(logging.DEBUG)
# log.addHandler(logging.StreamHandler())


class SnpeffEffect():
    hgvsparser = hgvs.parser.Parser()

    def __str__(self):
        return "Snpeff(%s; %s; %s; %s;)" % (
            self.geneid, str(self.annotation), str(self.hgvs_c), str(self.hgvs_p) if self.hgvs_p else "")

    def __repr__(self):
        return self.__str__()

    def __init__(self, alt, effects, impact, gene, geneid, feature_type, feature_id,
                 transcript_biotype, rank_div_total, hgvs_c, hgvs_p, c_dna_pos,
                 cds_pos, aa_pos, dist_to_feature, errors,
                 aa_len):
        self.alt = alt
        self.annotation = effects
        self.impact = impact
        self.gene = gene
        self.geneid = geneid
        self.feature_type = feature_type
        self.feature_id = feature_id
        self.transcript_biotype = transcript_biotype
        self.rank_div_total = rank_div_total
        self.hgvs_c = hgvs_c
        self.hgvs_p = hgvs_p
        self.c_dna_pos = c_dna_pos
        self.cds_pos = cds_pos
        self.aa_pos = aa_pos
        self.dist_to_feature = int(dist_to_feature) if dist_to_feature else None
        self.errors = errors
        self.aa_len = aa_len
        self.aa_ref = ""
        self.aa_alt = ""
        if self.hgvs_c:
            self.gene_pos = hgvs_c.pos.start.base
        if self.hgvs_p and self.hgvs_p.pos:

            self.aa_ref = self.hgvs_p.pos.start.aa
            try:
                if self.hgvs_p.edit.type == "del":
                    self.aa_alt = "del"
                elif self.hgvs_p.edit.type == "dup":
                    self.aa_alt = "dup"
                else:
                    self.aa_alt = self.hgvs_p.edit.alt
            except:
                _log.warn(self.hgvs_p.edit)

    @classmethod
    def read(cls, ann_str):

        aa_pos, aa_len = (None, None)
        (alt, annotation, impact, gene, geneid, feature_type, feature_id, transcript_biotype,
         rank_div_total, hgvs_c, hgvs_p, c_dna_pos, cds_pos, (aa_pos_aa_len), dist_to_feature, errors) = ann_str.split(
            "|")
        annotation = annotation.split("&")
        if aa_pos_aa_len:
            aa_pos, aa_len = aa_pos_aa_len.split("/")
            aa_pos = int(aa_pos)
        try:
            hgvs_c = cls.hgvsparser.parse_hgvs_variant("xx:" + hgvs_c).posedit
        except Exception as ex:
            _log.warn(ex)
            hgvs_c = ""
        if hgvs_p:
            try:
                hgvs_p = cls.hgvsparser.parse_hgvs_variant("xx:" + hgvs_p).posedit
            except Exception as ex:
                _log.warn(ex)
                hgvs_p = ""
        else:
            hgvs_p = None

        return SnpeffEffect(alt, annotation, impact, gene, geneid, feature_type, feature_id, transcript_biotype,
                            rank_div_total, hgvs_c, hgvs_p, c_dna_pos, cds_pos, aa_pos, dist_to_feature, errors, aa_len)


class VcfSnpeffIO():

    @classmethod
    def parse(cls, vcf_path):
        if hasattr(vcf_path, "read"):
            h = vcf_path
        else:
            if vcf_path.endswith(".gz"):
                h = gzip.open(vcf_path)
            else:
                h = open(vcf_path)

        try:
            variantes = vcf.VCFReader(h)
            for v in variantes:
                effects = [SnpeffEffect.read(x) for x in (v.INFO["ANN"] if "ANN" in v.INFO else [])]
                intergenic = [(i, x) for i, x in enumerate(effects) if "intragenic_variant" in x.annotation]
                if intergenic:
                    i, intergenic = intergenic[0]
                    if (("upstream_gene_variant" in effects[0].annotation)
                            or ("downstream_gene_variant" in effects[0].annotation)):
                        effects = [effects[i]] + effects[:i - 1] + effects[i:]
                yield (v, effects)
        finally:
            h.close()


class Command(BaseCommand):
    help = 'Loads the obo files to the database. Prepared for GO and SO'

    def __init__(self, stdout=None, stderr=None, no_color=False):
        super().__init__(stdout=stdout, stderr=stderr, no_color=no_color)
        self.gene_cache = {}
        self.contig_cache = {}

    def search_gene(self, ref, gid):
        if gid in self.gene_cache:
            gene = self.gene_cache[gid]
        else:
            gene = Seqfeature.objects.get(Q(bioentry=ref, type_term__name="gene") &
                                          Q(qualifiers__value=gid,
                                            qualifiers__term__name="locus_tag"))
            self.gene_cache[gid] = gene
        return gene

    # def search_gene_pos(self, ref, pos):
    #     gene = Seqfeature.objects.filter(bioentry=ref, type_term__name="gene", locations__start_pos__lte=pos,
    #                                      locations__end_pos__gt=pos)
    #     if gene.exists():
    #         gene = gene.get()
    #     else:
    #         gene_fw = Seqfeature.objects.filter(bioentry=ref, type_term__name="gene", locations__strand=1,
    #                                             locations__start_pos__gt=pos).order_by("locations__start_pos")
    #         gene_rv = Seqfeature.objects.filter(bioentry=ref, type_term__name="gene", locations__strand=-1,
    #                                             locations__start_pos__lt=pos).order_by("locations__start_pos")
    #         if gene_fw.exists() and gene_rv.exists():
    #             gene_fw = gene_fw[0]
    #             gene_rv = gene_rv[gene_rv.count()-1]
    #             gene = gene_fw if (abs(gene_fw.locations.all()[0].start_pos - pos) < abs(
    #                 (gene_rv.locations.all()[0].end_pos - pos))) else gene_rv
    #         elif gene_fw.exists():
    #             gene = gene_fw[0]
    #         else:
    #             gene = gene_rv[gene_rv.count()-1]
    #
    #     return gene

    def select_effect(self, effects, alt):
        allele_effects = [x for x in effects if x.alt == alt]
        allele_effects2 = allele_effects
        intergenic = [(i, x) for i, x in enumerate(allele_effects) if "intragenic_variant" in x.annotation]
        if intergenic:
            i, intergenic = intergenic[0]
            if (("upstream_gene_variant" in allele_effects[0].annotation)
                    or ("downstream_gene_variant" in allele_effects[0].annotation)):
                allele_effects2 = [allele_effects[i]] + allele_effects[:i - 1] + allele_effects[i:]
        if len(set(reduce(lambda x, y: x + y, [x.annotation for x in allele_effects])) - set(
                ["upstream_gene_variant", "downstream_gene_variant", "intergenic_region"])) == 0:

            allele_effects2 = sorted(
                [x for x in allele_effects if x.dist_to_feature and "upstream_gene_variant" in x.annotation],
                key=lambda x: x.dist_to_feature)
            if not allele_effects2:
                allele_effects2 = sorted([x for x in allele_effects if x.dist_to_feature],
                                        key=lambda x: x.dist_to_feature)
        return allele_effects2[0]

    def add_allele(self, sample, vc, variant, new_variant, effects):
        if sample.data.GT in ["CN" + x for x in [""] + [str(y) for y in range(20)]]:
            self.stderr.write("GT ERROR %s" % sample.data.GT)
            return

        if sample.called:
            alt = str(variant.ALT[int(sample.data.GT) - 1])
        elif sample.data.GT == ".":
            return
        else:
            raise Exception("Not expected...")

        effect = self.select_effect(effects, alt)

        allele_query = Allele.objects.prefetch_related("variant_fk").filter(variant_fk=new_variant, alt=alt)
        if not allele_query.exists():
            new_allele = Allele(variant_fk=new_variant, alt=alt)
            self.add_effect(effect, new_allele, new_variant, vc)
            new_allele.hgvs_c = str(effect.hgvs_c)
            new_allele.save()
        else:
            new_allele = allele_query.get()
            assert set(new_allele.main_effect_fk.variant_type.split("|")) & set(effect.annotation)

        # assert new_variant.gene_pos == effect.gene_pos, [new_variant.gene_pos, effect.gene_pos]
        # assert abs(abs(new_variant.gene_pos) - abs(effect.gene_pos))<= max(len(variant.REF),len(new_allele.alt)), [new_variant.gene_pos, effect.gene_pos]

        if not Variantassignment.objects.filter(variant_collection_fk=vc, variant_fk=new_variant,
                                                allele_fk=new_allele).exists():
            assignment = Variantassignment(variant_collection_fk=vc, variant_fk=new_variant, allele_fk=new_allele)
            assignment.save()

            assert sample.data.GT != "."
            for k, v in sample.data._asdict().items():
                val = "|".join(map(str, v)) if isinstance(v, (list, tuple)) else str(v)
                Variantannotation(
                    source_type="prediction", source="VC",
                    assignment_fk=assignment, prop=k, value=val).save()

            qual_va = Variantannotation(
                source_type="prediction", source="VC",
                assignment_fk=assignment, prop="qual", value=variant.QUAL)
            qual_va.save()

    def add_effect(self, effect, new_allele, new_variant, vc):

        gene = self.search_gene(new_variant.contig, effect.feature_id)
        if effect.aa_pos:
            effect_query = Effect.objects.filter(transcript=effect.geneid, ref_organism=vc.ref_organism,
                                                 variant_type="|".join(effect.annotation),
                                                 aa_pos=effect.aa_pos, aa_ref=effect.aa_ref, aa_alt=effect.aa_alt)
            if effect_query.exists():
                new_effect = effect_query.get()
            else:
                new_effect = Effect(transcript=effect.geneid, ref_organism=vc.ref_organism,
                                    variant_type="|".join(effect.annotation), gene=gene)
                new_effect.aa_pos = effect.aa_pos
                new_effect.aa_ref = effect.aa_ref
                new_effect.aa_alt = effect.aa_alt
                new_effect.hgvs_p = str(effect.hgvs_p)
        else:
            new_effect = Effect(transcript=effect.geneid, ref_organism=vc.ref_organism,
                                variant_type="|".join(effect.annotation), gene=gene)
        assert gene == new_variant.gene
        new_effect.save()
        new_allele.main_effect_fk = new_effect
        new_allele.save()

        if not AlleleEffect.objects.filter(allele_fk=new_allele, effect_fk=new_effect).exists():
            AlleleEffect(allele_fk=new_allele, effect_fk=new_effect).save()

        return new_effect

    def add_variant(self, ref_organism, samples_dict, var, effects):
        if var.CHROM in self.contig_cache:
            seq = self.contig_cache[var.CHROM]
        else:
            seq = Bioentry.objects.get(biodatabase=ref_organism, identifier=var.CHROM)
            self.contig_cache[var.CHROM] = seq
        variant_query = Variant.objects.filter(
            pos=var.POS, contig=seq, ref=var.REF)

        if not variant_query.exists():
            # gene = self.search_gene_pos(seq, var.POS)
            effect = self.select_effect(effects,str(var.ALT[int(var.samples[0].data.GT) - 1]))
            gene = self.search_gene(seq, effect.feature_id)

            genepos = (
                          var.POS - gene.locations.all()[0].start_pos + 1
                          if gene.locations.all()[0].strand == 1
                          else gene.locations.all()[0].end_pos - var.POS + 1)

            if genepos < 0:
                genepos = genepos - 1

            new_variant = Variant(contig=seq, pos=var.POS, gene_pos=genepos, ref=var.REF, gene=gene)
            new_variant.save()

        else:
            new_variant = variant_query.get()

        for sample in var.samples:
            self.add_allele(sample, samples_dict[sample.sample], var, new_variant, effects)

    def exists_sample(self, ref_organism, sample):
        return Variantcollection.objects.filter(
            sample=sample, ref_organism=ref_organism).exists()

    def delete_sample(self, ref_organism, sample):
        return Variantcollection.objects.where(
            sample=sample, ref_organism=ref_organism).delete()

    def load_variants(self, variants, ref_organism, samples, total):
        '''
        variants: array of (variant,effects ) tuple
        '''
        samples_dict = {}

        for sample in samples:
            sample_name = sample.sample.split(".variant")[0]
            if self.exists_sample(ref_organism, sample_name):
                # raise VariantcollectionExistsError(ref_organism, sample.sample)
                vc = Variantcollection.objects.get(ref_organism=ref_organism, sample=sample_name)
            else:
                vc = Variantcollection(ref_organism=ref_organism, sample=sample_name)
                vc.save()
            samples_dict[sample.sample] = vc

        process = []
        for idx, (var, effects) in enumerate(tqdm(variants, total=total), 1):
            process.append([ref_organism, samples_dict, var, effects])
            if len(process) == 100:
                with transaction.atomic():
                    for x in tqdm(process):
                        self.add_variant(*x)
                process = []
        with transaction.atomic():
            for i,x in tqdm(enumerate(process)):
                self.add_variant(*x)

    def add_arguments(self, parser):
        parser.add_argument('--vcf', required=True)
        parser.add_argument('--reference', required=True)

    def handle(self, *args, **options):
        # try:
        if options["vcf"].endswith(".gz"):
            h = gzip.open(options["vcf"], 'r')
        else:
            h = open(options["vcf"])
        variant = next(vcf.VCFReader(h))
        # finally:
        #     h.close()

        ref_organism = Biodatabase.objects.get(name=options["reference"])
        total = int(subprocess.check_output("grep -v ^#  " + options["vcf"] + " | wc -l", shell=True).split()[0])
        self.load_variants(VcfSnpeffIO.parse(options["vcf"]), ref_organism, variant.samples, total)
