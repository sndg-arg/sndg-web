from collections import defaultdict


class PhenoGenoTable():
    def __init__(self):
        self._sample_pheno = []
        self.stats = {}
        self.stats_genetic = {}
        self._pheno_sample = {}

    def sample_pheno(self, vca):
        self._sample_pheno = []
        reported_alleles = {}

        for vc in vca.samples():
            sample = {"name": vc.sample, "phenos": []}

            for p in vca.phenotypes():

                pheno = {}
                gquery = [g for g in vc.genotypes.all() if g.phenotype == p]

                if gquery:
                    g = gquery[0]
                    pheno["genotype"] = list(g.supportedby.all())
                    reported_alleles.update(
                        {x.reported_allele.allele: x.reported_allele for x in pheno["genotype"] if x.reported_allele})

                else:
                    pheno["genotype"] = []

                aquery = [a for a in vc.assays.all() if a.protocol.phenotype == p]
                if aquery:
                    a = aquery[0]
                    pheno["phenotype"] = a.result.identifier
                else:
                    pheno["phenotype"] = "?"
                sample["phenos"].append(pheno)

            self._sample_pheno.append(sample)
        return self._sample_pheno

    def pheno_sample(self, collection_set):
        """

        :param collection_set: Collection_set
        :return:
        Dict<pheno_name,{name,samples:[{"sample":..,genotype:[ReportedAllele,...]},...]} >
        """
        if self._pheno_sample:
            return self._pheno_sample

        self._pheno_sample = []
        reported_alleles = {}

        for p in collection_set.phenotypes():
            pheno = {"name": p.name, "samples": []}

            for vc in collection_set.samples():

                sample = {"sample": vc.sample}
                gquery = [g for g in vc.genotypes.all() if g.phenotype == p]

                if gquery:
                    g = gquery[0]
                    sample["genotype"] = list(g.supportedby.all())
                    reported_alleles.update(
                        {x.reported_allele.allele: x.reported_allele for x in sample["genotype"] if x.reported_allele})

                else:
                    sample["genotype"] = []
                aquery = [a for a in vc.assays.all() if a.protocol.phenotype == p]
                if aquery:
                    a = aquery[0]
                    sample["phenotype"] = a.result.identifier
                else:
                    sample["phenotype"] = "?"
                pheno["samples"].append(sample)

            self._pheno_sample.append(pheno)
        return self._pheno_sample

    def _get_strain(self, sample_name):
        return [x for x in self.variant_collection_set.samples() if x.sample == sample_name][0]

    def pheno_sample_stats(self, variant_collection_set, pos, neg, conclusive, possible, hint,
                           geno_status_filter=lambda g: True):
        self.variant_collection_set = variant_collection_set
        self._pheno_sample = self.pheno_sample(variant_collection_set)
        self.stats = {x["name"]: {} for x in self._pheno_sample}
        self.stats_genetic = {x["name"]: defaultdict(lambda: defaultdict(lambda: [])) for x in self._pheno_sample}

        gtpt_evaluations = len(variant_collection_set.studied_phenotypes.all()) * len(
            variant_collection_set.assignments.all())
        for p in self._pheno_sample:
            for s in p["samples"]:
                if s["phenotype"] == "?":
                    gtpt_evaluations -= 1

        for pheno_dict in self._pheno_sample:
            self._process_phenotype(conclusive, geno_status_filter, hint, neg, pheno_dict, pos, possible)

        self.stats["total"] = {
            x: sum([(len(self.stats[p][x]) if len(x) == 2 else self.stats[p][x]) for p in self.stats])
            for x in ["TP", "TN", "FP", "FN"]}
        with_phenotype = (self.stats["total"]["TP"] + self.stats["total"]["FN"])
        if with_phenotype:
            self.stats["total"]["sensitivity"] = self.stats["total"]["TP"] * 1.0 / with_phenotype
        else:
            self.stats["total"]["sensitivity"] = 0
        without_phenotype = (self.stats["total"]["TN"] + self.stats["total"]["FP"])
        if without_phenotype:
            self.stats["total"]["specificity"] = self.stats["total"]["TN"] * 1.0 / without_phenotype
        else:
            self.stats["total"]["specificity"] = 0
        self.stats["total"]["tested"] = gtpt_evaluations
        return self.stats, self.stats_genetic, self._pheno_sample

    def _process_phenotype(self, conclusive, geno_status_filter, hint, neg, pheno_dict, pos, possible):
        pheno = pheno_dict["name"]
        samples = pheno_dict["samples"]

        no_assays = [sample for sample in samples if len([x for x in self._get_strain(sample["sample"]).assays.all() if
                                                          x.protocol.phenotype.name == pheno_dict["name"]]) == 0]

        self.stats[pheno]["TP"] = {}
        for x in samples:
            reported_variants = [gs for gs in x["genotype"] if gs.status in [conclusive, possible]]
            if (x["phenotype"] == pos.identifier) and reported_variants:
                self.stats[pheno]["TP"][x["sample"]] = reported_variants
                for rv in reported_variants:
                    gene = rv.assignment.allele_fk.variant_fk.gene.locus_tag()
                    if geno_status_filter(rv) and (
                            x["sample"] not in self.stats_genetic[pheno][gene + " " + rv.assignment.allele_fk.hgvs_c][
                        "TP"]):
                        self.stats_genetic[pheno][gene + " " + rv.assignment.allele_fk.hgvs_c]["TP"].append(x["sample"])
        self.stats[pheno]["TN"] = {}
        for x in samples:
            hinted_variants = [gs for gs in x["genotype"] if gs.status == hint]
            if (x["phenotype"] == neg.identifier) and (
                    len([gs for gs in x["genotype"] if gs.status in [conclusive, possible]]) == 0):
                self.stats[pheno]["TN"][x["sample"]] = hinted_variants
                for rv in hinted_variants:
                    gene = rv.assignment.allele_fk.variant_fk.gene.locus_tag()
                    if geno_status_filter(rv) and (
                            x["sample"] not in self.stats_genetic[pheno][gene + " " + rv.assignment.allele_fk.hgvs_c][
                        "TN"]):
                        self.stats_genetic[pheno][gene + " " + rv.assignment.allele_fk.hgvs_c]["TN"].append(x["sample"])
        self.stats[pheno]["FP"] = {}
        for x in samples:
            reported_variants = [gs for gs in x["genotype"] if gs.status in [conclusive, possible]]
            if (x["phenotype"] == neg.identifier) and reported_variants:
                self.stats[pheno]["FP"][x["sample"]] = reported_variants
                for rv in reported_variants:
                    gene = rv.assignment.allele_fk.variant_fk.gene.locus_tag()
                    if geno_status_filter(rv) and (
                            x["sample"] not in self.stats_genetic[pheno][gene + " " + rv.assignment.allele_fk.hgvs_c][
                        "FP"]):
                        self.stats_genetic[pheno][gene + " " + rv.assignment.allele_fk.hgvs_c]["FP"].append(x["sample"])
        self.stats[pheno]["FN"] = {}
        for x in samples:
            hinted_variants = [gs for gs in x["genotype"] if gs.status == hint]
            if (x["phenotype"] == pos.identifier) and (
                    len([gs for gs in x["genotype"] if gs.status in [conclusive, possible]]) == 0):
                self.stats[pheno]["FN"][x["sample"]] = hinted_variants
                for rv in hinted_variants:
                    gene = rv.assignment.allele_fk.variant_fk.gene.locus_tag()
                    if geno_status_filter(rv) and (
                            x["sample"] not in self.stats_genetic[pheno][gene + " " + rv.assignment.allele_fk.hgvs_c][
                        "FN"]):
                        self.stats_genetic[pheno][gene + " " + rv.assignment.allele_fk.hgvs_c]["FN"].append(x["sample"])
        resistant = len(self.stats[pheno]["TP"]) + len(self.stats[pheno]["FN"])
        if resistant:
            self.stats[pheno]["sensitivity"] = len(self.stats[pheno]["TP"]) * 1.0 / resistant
        else:
            self.stats[pheno]["sensitivity"] = 0
        sensible = len(self.stats[pheno]["TN"]) + len(self.stats[pheno]["FP"])
        if sensible:
            self.stats[pheno]["specificity"] = len(self.stats[pheno]["TN"]) * 1.0 / sensible
        else:
            self.stats[pheno]["specificity"] = 0

        # Chequeo set([x["sample"] for x in samples]) - set( reduce(    for x in ["TP","TN","FN","FP"], lambda x,y:x+y) )
        tot = []
        for x in ["TP", "TN", "FN", "FP"]:
            tot += self.stats[pheno][x].keys()
        assert len(tot) == len(set(tot))

        assert (len(tot) + len(no_assays)) == len(samples), [len(tot), len(no_assays), len(samples)]
