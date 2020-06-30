import os
import warnings
from collections import defaultdict

from Bio import BiopythonWarning, BiopythonParserWarning, BiopythonDeprecationWarning, BiopythonExperimentalWarning
from django.core.management.base import BaseCommand, CommandError
from tqdm import tqdm

warnings.simplefilter('ignore', RuntimeWarning)
warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)
warnings.simplefilter('ignore', BiopythonDeprecationWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

from bioseq.models.Ontology import Ontology

from SNDG.Sequence import smart_parse

"""
https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/
https://github.com/ncbi/pgap
https://www.ncbi.nlm.nih.gov/genbank/genome_validation/
https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/reannotation/
https://www.ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/
https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/
https://www.ncbi.nlm.nih.gov/genome/annotation_prok/
https://www.ebi.ac.uk/ena/submit/annotation-checklists
"""


class _SeqCheck:

    def process_entry(self, seqrecord):
        pass

    def process_feature(self, seqrecord, seqfeature, parent_seqfeature, children):
        pass

    def summarize(self):
        pass


class CheckRepeatedOrEmptyContig(_SeqCheck):

    def __init__(self):
        self.features_per_contig = defaultdict(lambda: 0)
        self.errors = []
        self.contig_names = {}
        self.repeated = False

    def process_entry(self, seqrecord):
        if len(seqrecord.seq) < 100:
            self.errors.append(f"{seqrecord.id} short contig {len(seqrecord.seq)}")
        if seqrecord.id in self.contig_names:
            self.errors.append(f"{seqrecord.id} is repeated")
            self.repeated = True
        self.contig_names[seqrecord.id] = 1
        self.features_per_contig[seqrecord.id] = len(seqrecord.features)

    def summarize(self):
        pass


class CheckRepeatedOrNoLocusTag(_SeqCheck):
    def __init__(self):
        self.lts = {"genes": [], "other": []}
        self.repeated_lt = False
        self.errors = []
        self.errors_ex = {}
        self.errors_count = defaultdict(lambda: 0)

    def process_feature(self, seqrecord, seqfeature, parent_seqfeature, children):
        if seqfeature.type in ["source", "region","regulatory","repeat_region"]:
            return

        if "locus_tag" not in seqfeature.qualifiers:
            self.errors_ex = seqfeature.qualifiers
            self.errors_count["no locus_tag"] += 1
        else:
            if seqfeature.type == Ontology.gene:
                lt = seqfeature.qualifiers["locus_tag"][0]
                if lt in self.lts:
                    self.errors.append(f"{lt} repeated in different genes")
                    self.errors_count["lt_repeated"] += 1
                    self.errors_ex["lt_repeated"] += seqfeature

            else:
                lt = seqfeature.qualifiers["locus_tag"][0]
                if lt in self.lts["other"]:
                    if lt in self.lts["other"]:
                        self.errors.append(f"{lt} repeated in different products")
                        self.errors_count["lt_repeated"] += 1
                        self.errors_ex["lt_repeated"] += seqfeature
                    self.lts["other"][lt] = 1


class CheckInvalidTypes(_SeqCheck):

    def __init__(self):
        self.errors_count = defaultdict(lambda: 0)
        self.errors_ex = {}

    def process_feature(self, seqrecord, seqfeature, parent_seqfeature, children):
        if seqfeature.type == Ontology.gene:
            pass
        if seqfeature.type == Ontology.mRNA:
            self.gene_parent(parent_seqfeature)
            if len(children) == 0:
                err = f"{Ontology.mRNA} has no CDS associated"
                self.errors_count[err] += 1
                self.errors_ex[err] = seqfeature

        if seqfeature.type in [Ontology.rRNA, Ontology.tRNA, Ontology.misc_RNA,
                               Ontology.miRNA, Ontology.ncRNA]:
            self.gene_parent(parent_seqfeature)

        if seqfeature.type == Ontology.CDS:
            self.gene_parent(parent_seqfeature)

        if seqfeature.type == Ontology.mat_peptide:
            self.gene_parent(parent_seqfeature)

    def gene_parent(self, parent_seqfeature):
        if parent_seqfeature.type != "gene":
            self.errors_count[f"{Ontology.mRNA} parent is not a gene"] += 1


class CheckOverlapps(_SeqCheck):
    pass


class CheckCDSTranslation(_SeqCheck):
    def __init__(self):
        self.errors = defaultdict(lambda: [])
        self.errors_ex = {}
        self.errors_count = defaultdict(lambda: 0)

    def process_feature(self, seqrecord, seqfeature, parent_seqfeature, children):
        if seqfeature.type == Ontology.CDS:
            if "translation" not in seqfeature.qualifiers:
                if "pseudo" not in seqfeature.qualifiers:
                    self.errors_count["no translation"] += 1
            else:
                translated = seqfeature.extract(seqrecord.seq).translate()
                if str(translated) == seqfeature.qualifiers["translation"][0]:
                    self.errors["translation error"].append(seqfeature)


class CheckECs(_SeqCheck):
    pass


class Command(BaseCommand):
    help = 'checks if genebank meets the requirements to be loaded in the db '

    def add_arguments(self, parser):
        parser.add_argument('--input', '-i', required=True)
        parser.add_argument('--fix', nargs='+')

    def handle(self, *args, **options):
        if not os.path.exists(options['input']):
            raise CommandError(f'{options["input"]} does not exists')

        self.stderr.write(f"cheking  {options['input']} imported!")

        checks = [CheckRepeatedOrEmptyContig(), CheckInvalidTypes(), CheckCDSTranslation(),
                  CheckRepeatedOrNoLocusTag()]

        with tqdm(smart_parse(options["input"]), file=self.stderr, position=0, leave=True) as pbar:
            for contig in pbar:
                for check in checks:
                    check.process_entry(contig)

                    for i, f in enumerate(contig.features):
                        parent = contig.features[i - 1] if i > 0 else None
                        children = contig.features[i + 1] if i < (len(contig.features) - 1) else None
                        check.process_feature(contig, f, parent, children)

        for check in checks:
            print(check.__dict__)

        self.stderr.write(f"genome {options['input']} imported!")
