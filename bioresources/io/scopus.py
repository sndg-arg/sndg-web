import json

from elsapy.elsclient import ElsClient
from elsapy.elssearch import ElsSearch
from django.conf import settings


class ScopusDS:
    """
    client.inst_token = config['insttoken']
    https://dev.elsevier.com/documentation/ScopusSearchAPI.wadl#d1e33
    https://dev.elsevier.com/guides/ScopusSearchViews.htm
    """

    words = ["genome", "genomic", "genomics", "assembly", "assemblies", "transcriptome", "transcriptomics",
             "proteome", "proteomics", "sequencing", "resistome",
             "metabolome", "kinome", "reactome",
             "wgs", "rnaseq", "rna-seq", "drugome", "variome", "secretome", "pangenome", "metagenomic", "metagenomics",
             "differential expression", "gene prediction", "gene annotation",
             "regulome", "regulation network", "transposon", "transposase", "repeated regions", "scaffold", "contig",
             "assembler",
             "reads mapping", "operon", "genic island", "scaffold", "contig", "gene enhancer", "regulation motif",
             "exon",
             "intron",
             "alternative splicing", "illumina", "expression array", "chip-seq", "gwas", "snps"
             ]

    def __init__(self, con_path):
        with open(con_path) as con_file:
            self.config = json.load(con_file)

        self.client = ElsClient(self.config['apikey'])

    def query(self, country, keywords=words, after_year=None):
        """

        Example:
        aff_srch = ElsSearch('( ( AFFILCOUNTRY ( argentina )  AND  (TITLE-ABS-KEY ( Puccinia ) OR TITLE-ABS-KEY(Lactobacillus)) ) )  AND  ( burguener ) ','scopus')
        aff_srch.execute(client)
         print ("aff_srch has", len(aff_srch.results), "results.")
         AFFILCOUNTRY ( argentina )  AND  ( LIMIT-TO ( SRCTYPE ,  "j" )  OR  LIMIT-TO ( SRCTYPE ,  "p" ) )
        :return:
        """

        kwfilter = "OR ".join(['TITLE-ABS-KEY ("' + w + '") ' for w in keywords])

        after_year_q = "AND PUBYEAR AFT " + str(after_year) if after_year else ""

        aff_srch = ElsSearch("((" + kwfilter + ") AND AFFILCOUNTRY ( " + country + " ) " + after_year_q + " )",
                             'scopus', maxResults=10000)
        aff_srch.execute(self.client, get_all=True)

        for article in aff_srch.results:
            yield article

    def doi(self, doi):

        aff_srch = ElsSearch("DOI ( %s )" % doi, 'scopus', maxResults=1)
        aff_srch.execute(self.client, get_all=True)

        rs = [x for x in aff_srch.results if "error" not in x]
        if rs:
            doi = rs[0]['prism:doi']
            return {"doi": doi, "title": rs[0]["dc:title"], "record": rs[0]}
        else:
            return None


def scopus_client():
    return ScopusDS(settings.SCOPUS_API)
