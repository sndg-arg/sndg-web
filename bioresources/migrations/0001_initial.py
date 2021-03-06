# Generated by Django 2.2.2 on 2019-08-24 13:59

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('bioseq', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Organization',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=200)),
                ('url', models.URLField(null=True)),
                ('country', models.CharField(max_length=200, null=True)),
                ('city', models.CharField(max_length=200, null=True)),
                ('scopus_id', models.CharField(max_length=200, null=True)),
                ('scopus_names', models.TextField(null=True)),
                ('deprecated', models.BooleanField(default=False)),
                ('index_updated', models.BooleanField(default=False)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
            ],
            options={
                'verbose_name_plural': 'Organizations',
            },
        ),
        migrations.CreateModel(
            name='Person',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('surname', models.CharField(max_length=200)),
                ('name', models.CharField(default='', max_length=200)),
                ('scopus_id', models.CharField(max_length=200, null=True)),
                ('scopus_names', models.TextField(null=True)),
                ('email', models.EmailField(max_length=254)),
                ('deprecated', models.BooleanField(default=False)),
                ('index_updated', models.BooleanField(default=False)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
            ],
            options={
                'verbose_name_plural': 'Persons',
            },
        ),
        migrations.CreateModel(
            name='Resource',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('type', models.PositiveIntegerField(choices=[(0, 'PUBLICATION'), (1, 'BIOPROJECT'), (2, 'SEQUENCE'), (3, 'ASSEMBLY'), (4, 'GENOME'), (5, 'READS'), (6, 'STRUCTURE'), (7, 'EXPRESSION'), (8, 'BARCODE'), (9, 'SAMPLE'), (10, 'TOOL'), (40, 'PROTEIN'), (30, 'ORGANIZATION'), (20, 'PERSON')])),
                ('name', models.CharField(max_length=350)),
                ('description', models.TextField(blank=True)),
                ('deprecated', models.BooleanField(default=False)),
                ('index_updated', models.BooleanField(default=False)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('collaborators', models.ManyToManyField(blank=True, related_name='resources', to='bioresources.Person')),
                ('creators', models.ManyToManyField(blank=True, related_name='created_resources', to='bioresources.Organization')),
            ],
        ),
        migrations.CreateModel(
            name='ResourceProperty',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('organization', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioresources.Organization')),
                ('resource', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='properties', to='bioresources.Resource')),
                ('term', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Term')),
            ],
        ),
        migrations.CreateModel(
            name='RKeyword',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=200, unique=True)),
            ],
        ),
        migrations.CreateModel(
            name='Assembly',
            fields=[
                ('resource_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='bioresources.Resource')),
                ('intraspecific_name', models.CharField(max_length=250, null=True)),
                ('species_name', models.CharField(max_length=200, null=True)),
                ('level', models.PositiveIntegerField(choices=[(0, 'complete genome'), (1, 'chromosome'), (2, 'scaffold'), (3, 'contig')], null=True)),
                ('ncbi_org', models.CharField(max_length=200, null=True)),
                ('release_date', models.DateField(null=True)),
                ('update_date', models.DateField(null=True)),
                ('assembly_type', models.PositiveIntegerField(choices=[(0, 'haploid'), (1, 'diploid'), (2, 'other'), (3, 'unresolved-diploid')], null=True)),
            ],
            options={
                'verbose_name_plural': 'Assemblies',
            },
            bases=('bioresources.resource',),
        ),
        migrations.CreateModel(
            name='Barcode',
            fields=[
                ('resource_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='bioresources.Resource')),
                ('country', models.CharField(max_length=100)),
                ('subdivision', models.CharField(max_length=150)),
                ('marker', models.CharField(max_length=50, null=True)),
                ('image_url', models.URLField(null=True)),
                ('bold_org', models.CharField(max_length=255, null=True)),
                ('collectors', models.CharField(max_length=255, null=True)),
            ],
            options={
                'verbose_name_plural': 'Barcodes',
            },
            bases=('bioresources.resource',),
        ),
        migrations.CreateModel(
            name='BioProject',
            fields=[
                ('resource_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='bioresources.Resource')),
                ('sample_scope', models.PositiveIntegerField(choices=[(0, 'monoisolate'), (1, 'multi-species'), (2, 'environment'), (3, 'synthetic'), (4, 'other')], null=True)),
                ('material', models.PositiveIntegerField(choices=[(0, 'genome'), (1, 'metagenome'), (2, 'chromosome'), (3, 'transcriptome'), (4, 'reagent'), (5, 'proteome')], null=True)),
                ('capture', models.PositiveIntegerField(choices=[(0, 'whole'), (1, 'exome'), (2, 'barcode'), (3, 'TargetedLocusLoci')], null=True)),
                ('target', models.CharField(max_length=200, null=True)),
                ('submitters', models.TextField(null=True)),
                ('method', models.CharField(max_length=200, null=True)),
            ],
            options={
                'verbose_name_plural': 'BioProjects',
            },
            bases=('bioresources.resource',),
        ),
        migrations.CreateModel(
            name='Expression',
            fields=[
                ('resource_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='bioresources.Resource')),
                ('pdat', models.DateField(null=True)),
                ('gdstype', models.CharField(max_length=250, null=True)),
                ('submitters', models.TextField(null=True)),
            ],
            options={
                'verbose_name_plural': 'Expression',
            },
            bases=('bioresources.resource',),
        ),
        migrations.CreateModel(
            name='Publication',
            fields=[
                ('resource_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='bioresources.Resource')),
                ('doi', models.CharField(max_length=100)),
                ('date_of_publication', models.DateField(max_length=200)),
                ('pubmed_id', models.CharField(max_length=50, null=True)),
                ('electronic_id', models.CharField(max_length=50, null=True)),
                ('scopus_id', models.CharField(max_length=50, null=True)),
                ('issn', models.CharField(max_length=50, null=True)),
            ],
            options={
                'verbose_name_plural': 'Publications',
            },
            bases=('bioresources.resource',),
        ),
        migrations.CreateModel(
            name='ReadsArchive',
            fields=[
                ('resource_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='bioresources.Resource')),
                ('release_date', models.DateField(null=True)),
                ('update_date', models.DateField(null=True)),
            ],
            options={
                'verbose_name_plural': 'Reads Archive',
            },
            bases=('bioresources.resource',),
        ),
        migrations.CreateModel(
            name='Sample',
            fields=[
                ('resource_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='bioresources.Resource')),
                ('country', models.CharField(max_length=100)),
                ('subdivision', models.CharField(max_length=150)),
                ('collection_date', models.DateField(null=True)),
                ('publication_date', models.DateField(null=True)),
                ('update_date', models.DateField(null=True)),
            ],
            options={
                'verbose_name_plural': 'Samples',
            },
            bases=('bioresources.resource',),
        ),
        migrations.CreateModel(
            name='Structure',
            fields=[
                ('resource_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='bioresources.Resource')),
                ('pdbClass', models.CharField(max_length=50, null=True)),
                ('deposit_date', models.DateField(null=True)),
                ('method', models.CharField(max_length=50, null=True)),
                ('org_list', models.TextField(null=True)),
            ],
            options={
                'verbose_name_plural': 'Structures',
            },
            bases=('bioresources.resource',),
        ),
        migrations.CreateModel(
            name='Tool',
            fields=[
                ('resource_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='bioresources.Resource')),
                ('tool_type', models.PositiveSmallIntegerField(choices=[(0, 'app'), (1, 'database'), (2, 'library'), (3, 'plugin'), (4, 'program'), (5, 'webserver')])),
                ('url', models.URLField(null=True)),
            ],
            options={
                'verbose_name_plural': 'Tools',
            },
            bases=('bioresources.resource',),
        ),
        migrations.CreateModel(
            name='ResourcePropertyValue',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('value', models.CharField(max_length=200)),
                ('property', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, related_name='value', to='bioresources.ResourceProperty')),
            ],
        ),
        migrations.AddField(
            model_name='resource',
            name='keywords',
            field=models.ManyToManyField(related_name='associated_resources', to='bioresources.RKeyword'),
        ),
        migrations.AddField(
            model_name='resource',
            name='ncbi_tax',
            field=models.ForeignKey(db_column='ncbi_tax', null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='bioresources', to='bioseq.Taxon', to_field='ncbi_taxon_id'),
        ),
        migrations.AddField(
            model_name='resource',
            name='publishers',
            field=models.ManyToManyField(blank=True, related_name='published_resources', to='bioresources.Organization'),
        ),
        migrations.CreateModel(
            name='Job',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('command', models.TextField()),
                ('start', models.DateTimeField(null=True)),
                ('end', models.DateTimeField(null=True)),
                ('result', models.TextField(null=True)),
                ('status', models.PositiveIntegerField(choices=[(0, 'NEW'), (1, 'QUEUED'), (2, 'ERROR'), (3, 'RETRYING'), (4, 'FINISHED')], default=0)),
                ('retry', models.PositiveIntegerField(default=0)),
                ('dev_error', models.TextField(null=True)),
                ('user', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='jobs', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name_plural': 'Job',
            },
        ),
        migrations.CreateModel(
            name='Identity',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('identifier', models.CharField(max_length=200)),
                ('email', models.EmailField(max_length=254, null=True)),
                ('url', models.URLField(null=True)),
                ('authority', models.CharField(max_length=200)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField()),
                ('ends', models.DateTimeField()),
                ('person', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='bioresources.Person')),
            ],
            options={
                'verbose_name_plural': 'Identities',
            },
        ),
        migrations.CreateModel(
            name='ExternalId',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('identifier', models.CharField(max_length=20)),
                ('type', models.CharField(max_length=20)),
                ('organization', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='bioresources.Organization')),
                ('resource', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='external_ids', to='bioresources.Resource')),
            ],
            options={
                'verbose_name_plural': 'ExternalId',
            },
        ),
        migrations.CreateModel(
            name='Affiliation',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('author', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='affiliations', to='bioresources.Person')),
                ('organizations', models.ManyToManyField(related_name='affiliations', to='bioresources.Organization')),
                ('resource', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='affiliations', to='bioresources.Resource')),
            ],
            options={
                'verbose_name_plural': 'Affiliations',
            },
        ),
        migrations.CreateModel(
            name='ResourceRelation',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('role', models.CharField(max_length=200)),
                ('deprecated', models.BooleanField(default=False)),
                ('source', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='targets', to='bioresources.Resource')),
                ('target', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='sources', to='bioresources.Resource')),
            ],
            options={
                'verbose_name_plural': 'Resource Relations',
                'unique_together': {('source', 'target', 'deprecated')},
            },
        ),
        migrations.AlterUniqueTogether(
            name='resource',
            unique_together={('type', 'name')},
        ),
    ]
