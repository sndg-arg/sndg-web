import json
import sys

import Bio.SeqIO as bpio
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from chembl_model.models import TargetDictionary, ZincCompound, ZincProperty, Assays
from django.core.management.base import BaseCommand
from tqdm import tqdm


class Command(BaseCommand):
    help = 'Loads gene ontology terms'

    def add_arguments(self, parser):
        parser.add_argument('command')

    def target2json(self, protein_accession="P47791"):
        def tofloat(x):
            try:
                return float(x)
            except:
                return x
        data = {"target": {"accession": protein_accession, "chemblid": ""}, "assays": []}
        for a in Assays.objects.using("chembl_25").prefetch_related("activities__properties", "assay_type",
                                                                    "activities__record__molregno__chembl",
                                                                    "chembl","activities__record__molregno__properties").filter(
            tid__components__component__accession=protein_accession):

            data["target"]["chemblid"] = a.tid.chembl.chembl_id
            assay_data = {
                "chemblid": a.chembl.chembl_id,
                "type": a.assay_type.assay_desc,

                "description": a.description,
                "activities": []
            }
            data["assays"].append(assay_data)

            for act in a.activities.all():
                props = {
                    k:tofloat(v) for k,v in act.record.molregno.properties.__dict__.items()
                    if k != "_state"
                } if  hasattr(act.record.molregno,"properties") else {}
                act_data = {
                    "compound": {"key": act.record.compound_key, "name": act.record.compound_name,
                                 "chemblid": act.record.molregno.chembl.chembl_id,
                                 "is_drug": bool(act.record.drugs.all()),
                                 "properties":props},
                    "standard_relation": act.standard_relation,
                    "standard_value": tofloat( act.standard_value),
                    "standard_units": act.standard_units,
                    "activity_comment": act.activity_comment,
                    "pchembl_value": tofloat(act.pchembl_value),
                    "properties": [
                        {"type": x.type, "relation": x.relation, "value": float(x.value) if x.value else None,
                         "units": x.units, "text": x.text_value, "comments": x.comments} for x in act.properties.all()]
                }
                assay_data["activities"].append(act_data)
        return data

    def handle(self, *args, **options):
        if options["command"] == "fasta":
            """SELECT t.chembl_id AS target_chembl_id,
            t.pref_name        AS target_name,
            t.target_type,
            c.accession        AS protein_accession,
            c.sequence         AS protein_sequence
            FROM target_dictionary t
              JOIN target_type tt ON t.target_type = tt.target_type
              JOIN target_components tc ON t.tid = tc.tid
              JOIN component_sequences c ON tc.component_id = c.component_id
            AND tt.parent_type  = 'PROTEIN';
            """
            accessions = {}
            with sys.stdout as h:
                qs = TargetDictionary.objects.using("chembl_25").prefetch_related("components__component").filter(
                    target_type__parent_type="PROTEIN")
                total = qs.count()
                for target_dictionary in tqdm(qs, total=total, file=sys.stderr):
                    for target_component in target_dictionary.components.all():
                        rid = target_component.component.accession
                        if rid not in accessions:
                            accessions[rid] = 1
                            r = SeqRecord(id=rid, name="", description=target_dictionary.pref_name,
                                          seq=Seq(target_component.component.sequence))
                            bpio.write(r, h, "fasta")
        if options["command"] == "load_zinc":
            ZincCompound.objects.all().delete()
            ZincProperty.objects.all().delete()
            df = pd.read_table("/data/databases/zinc/6_prop.xls")
            bulk = []
            for i, r in tqdm(df.iterrows()):
                d = r.to_dict()
                d["code"] = d["ZINC_ID"]
                del d["ZINC_ID"]
                bulk.append(Zinc(**d))

                if i and (i % 1000 == 0):
                    Zinc.objects.bulk_create(bulk)
                    bulk = []
        if options["command"] == "dump_tp_assays":

            with open("/data/databases/target/processed/chembl_targets.txt") as h:
                targets = [x.strip() for  x in h]

            for x in tqdm(targets,file=sys.stderr):
                self.stdout.write(json.dumps(self.target2json(x)) + "\n")

