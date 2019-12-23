from django_tqdm import BaseCommand
from django.db.models import Q
from django.db import transaction

import os
import logging
from vardb.models import Allele, Variantannotation, Variantcollection, AntibioticResistance, VariantCollectionSet, Genotype, GenotypeSupport
import vcf
from tqdm import tqdm
import hgvs.parser

_log = logging.getLogger(__name__)


# log = logging.getLogger('django.db.backends')
# log.setLevel(logging.DEBUG)
# log.addHandler(logging.StreamHandler())

class Command(BaseCommand):
    help = 'Updates a variant collection genotypes'

    def __init__(self, stdout=None, stderr=None, no_color=False):
        super().__init__(stdout=stdout, stderr=stderr, no_color=no_color)

    def add_arguments(self, parser):
        parser.add_argument('-vc', '--variant_collection', help="variant collection sample name")
        parser.add_argument('-vs', '--collection_set', help="variant collection set name")

    def process_variant_collection(self, vc):
        pbar = tqdm(list(AntibioticResistance.objects.all()))
        for ar in pbar:
            pbar.set_description(str(ar))
            genotype = Genotype.objects.filter(phenotype=ar, variant_collection=vc)
            if genotype.exists():
                genotype = genotype.get()
            else:
                genotype = Genotype(phenotype=ar, variant_collection=vc)
                genotype.save()
            support = ar.process_variant_collection(genotype, vc)

            with transaction.atomic():
                for x in support:
                    if not GenotypeSupport.objects.filter(genotype=genotype,assignment=x.assignment).exists():
                        x.genotype = genotype
                        x.save()

    def handle(self, *args, **options):
        assert options["variant_collection"] or options[
            "collection_set"], "variant_collection or collection_set parameter must be specified"
        if options["variant_collection"]:
            self.process_variant_collection(Variantcollection.objects.get(sample=options["variant_collection"]))
        else:
            with tqdm(VariantCollectionSet.objects.get(name=options["collection_set"]).assignments.all()) as pbar:
                for vca in pbar:
                    pbar.set_description(vca.variant_collection.sample)
                    self.process_variant_collection(vca.variant_collection)
