from django_tqdm import BaseCommand

from biosql.models import Term, Ontology
from django_tqdm import BaseCommand

from vardb.models import Protocol, AntibioticResistance, GenotypeSupport, Assay


class Command(BaseCommand):
    help = 'Loads the obo files to the database. Prepared for GO and SO'

    def __init__(self, stdout=None, stderr=None, no_color=False):
        super().__init__(stdout=stdout, stderr=stderr, no_color=no_color)
        self.resist = None

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        """https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram-myco/
        https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/
         https://www.ncbi.nlm.nih.gov/biosample/docs/beta-lactamase/
        """
        # Ontology Assay Result
        if not Ontology.objects.filter(name=Assay.ASSAY_ONTOLOGY).exists():
            aro = Ontology(name=Assay.ASSAY_ONTOLOGY)
            aro.save()
            Term(name="Positive", identifier="Positive", ontology=aro).save()
            Term(name="Negative", identifier="Negative", ontology=aro).save()
            Term(name="Inconclusive", identifier="Inconclusive", ontology=aro).save()

        # Ontology Genotype Support Status
        if not Ontology.objects.filter(name=GenotypeSupport.STATUS_ONTOLOGY).exists():
            gsso = Ontology(name=GenotypeSupport.STATUS_ONTOLOGY)
            gsso.save()
            Term(name="Conclusive", identifier="Conclusive", ontology=gsso).save()
            Term(name="Possible", identifier="Possible", ontology=gsso).save()
            Term(name="Hint", identifier="Hint", ontology=gsso).save()

        # Ontology Antibiotic
        # Protocol : unknown...
        if not Ontology.objects.filter(name="Antibiotics").exists():
            ao = Ontology(name="Antibiotics")
            ao.save()
        else:
            ao = Ontology.objects.get(name="Antibiotics")

        for l in """streptomycin\tSTR
                    isoniazid\tINH
                    rifampicin\tRIF
                    ethambutol\tEMB
                    kanamycin\tKAN
                    amikacin\tAMK
                    capreomycin\tCAP
                    ethionamide\tETH
                    fluoroquinolone\tFLQ
                    aminoglycoside\tAMI
                    clofazimine\tCLO
                    bedaquilin\t
                    viomycin\tVIO
                    fosfomycin\tFOF
                    linezolid\tZD
                    pyrazinamide\tPZA                  
                    para-aminosalicylic_acid\tPAS""".split("\n"):
            name, ident = [x.strip() for x in l.split("\t")][:2]
            tquery = Term.objects.filter(name=name, identifier=ident, ontology=ao)
            if not tquery.exists():
                t = Term(name=name, identifier=ident, ontology=ao)
                t.save()
            else:
                t = tquery.get()
            aquery = AntibioticResistance.objects.filter(name="%s resistant" % name, antibiotic=t)
            if not aquery.exists():
                p = AntibioticResistance(name="%s resistant" % name, antibiotic=t)
                p.save()
            else:
                p = aquery.get()
            if not Protocol.objects.filter(name="%s antibiogram" % name, phenotype=p).exists():
                Protocol(name="%s antibiogram" % name, phenotype=p).save()








