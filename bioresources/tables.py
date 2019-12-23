import django_tables2 as tables

from .models import Publication


class PublicationTable(tables.Table):
    # counter = tables.TemplateColumn("{{ row_counter }}")

    class Meta:
        model = Publication
        template_name = 'django_tables2/bootstrap.html'
        exclude = ["type", "pubmed_id", "electronic_id", "scopus_id", "issn",
                   "ncbi_tax", "deprecated", "created_at", "updated_at","resource_ptr"]
