from neomodel import StructuredNode, StringProperty, DateProperty, One, RelationshipTo, IntegerProperty

from django_neomodel import DjangoNode
from neomodel import db

"""MATCH (a:Artist),(b:Album)
WHERE a.Name = "Strapping Young Lad" AND b.Name = "Heavy as a Really Heavy Thing"
CREATE (a)-[r:RELEASED]->(b)
RETURN r"""


def get_resource_class():
    from bioresources.models.Resource import Resource
    return {
        Resource.RESOURCE_TYPES.STRUCTURE: "Structure",
        Resource.RESOURCE_TYPES.ASSEMBLY: "Assembly",
        Resource.RESOURCE_TYPES.SAMPLE: "Sample",
        Resource.RESOURCE_TYPES.EXPRESSION: "Expression",
        Resource.RESOURCE_TYPES.PUBLICATION: "Publication",
        Resource.RESOURCE_TYPES.TOOL: "Tool",
        Resource.RESOURCE_TYPES.READS: "Reads",
        Resource.RESOURCE_TYPES.PERSON: "Person",
        Resource.RESOURCE_TYPES.ORGANIZATION: "Organization",
        Resource.RESOURCE_TYPES.BIOPROJECT: "BioProject",
        Resource.RESOURCE_TYPES.BARCODE: "Barcodes",

    }


def delete_edge(r1, r2, reltype: str = "USES"):
    resource_class = get_resource_class()
    query = """ MATCH (a:%(src)s)-[e:%(reltype)s]-(b:%(dst)s)
                WHERE a.rid = %(src_id)s AND b.rid = %(dst_id)s                 
                DELETE e""" % {
        "src": resource_class[r1.type],
        "dst": resource_class[r2.type],
        "src_id": r1.id, "dst_id": r2.id, "reltype": reltype
    }
    return db.cypher_query(query, {})


def connect_nodes(r1, r2, reltype: str = "USES"):
    resource_class = get_resource_class()

    query = """ MATCH (a:%(src)s),(b:%(dst)s)
                WHERE a.rid = %(src_id)s AND b.rid = %(dst_id)s 
                CREATE (a)-[r:%(reltype)s]->(b)
                RETURN r""" % {
        "src": resource_class[r1.type],
        "dst": resource_class[r2.type],
        "src_id": r1.id, "dst_id": r2.id, "reltype": reltype
    }
    return db.cypher_query(query, {})


class Tax(StructuredNode):
    rid = IntegerProperty(unique_index=True)
    name = StringProperty()
    ncbi_id = StringProperty(unique_index=True)


class Country(StructuredNode):
    name = StringProperty(unique_index=True)


class BioProject(StructuredNode):
    rid = IntegerProperty(unique_index=True)
    name = StringProperty(unique_index=True)

    # level = DateProperty()
    # location = RelationshipTo('Country', 'LOCATION')

    @classmethod
    def from_resource(cls, obj):
        n = cls(rid=obj.id, name=obj.name)
        n.save()
        return n


class Organization(StructuredNode):
    rid = IntegerProperty(unique_index=True)
    name = StringProperty(unique_index=True)
    level = DateProperty()
    location = RelationshipTo('Country', 'LOCATION')

    @classmethod
    def from_resource(cls, organization):
        o = cls(rid=organization.id, name=organization.name)
        o.save()
        return o


class Journal(StructuredNode):
    name = StringProperty(unique_index=True)


class Person(StructuredNode):
    rid = IntegerProperty(unique_index=True)
    name = StringProperty()

    @classmethod
    def from_resource(cls, person):
        gPerson = Person(rid=person.id, name=person.complete_name())
        gPerson.save()


class Species(StructuredNode):
    name = StringProperty(unique_index=True)


class Resource(StructuredNode):
    rid = IntegerProperty(unique_index=True)
    title = StringProperty(unique_index=True)
    # description = StringProperty()
    species = RelationshipTo('Species', 'SPECIES')
    tax = RelationshipTo('Tax', 'TAX')


class Publication(Resource):
    title = StringProperty(unique_index=True)
    # published = DateProperty()
    authors = RelationshipTo('Person', 'AUTHOR')
    organizations = RelationshipTo('Organization', 'AFF')
    journal = RelationshipTo('Journal', 'PUBLISHER', cardinality=One)
    resources = RelationshipTo('Resource', 'USED')

    @classmethod
    def from_resource(cls, publication):

        gPublication = Publication(rid=publication.id, title=publication.name,
                                   publication=publication.date_of_publication)
        gPublication.save()
        for aff in publication.affiliations.all():
            per = Person.nodes.get(rid=aff.author_id)
            gPublication.authors.connect(per)
            for org in aff.organizations.all():
                gOrg = Organization.nodes.get(rid=org.id)
                gPublication.organizations.connect(gOrg)
        return gPublication


class Expression(Resource):
    pdat = StringProperty()
    gdstype = StringProperty()

    @classmethod
    def from_resource(cls, expresion):
        r = Expression(rid=expresion.id, title=expresion.name, pdat=expresion.pdat, gdstype=expresion.gdstype)
        r.save()
        return r


class Assembly(Resource):
    intraspecific_name = StringProperty()
    level = StringProperty()
    reads = RelationshipTo('Reads', 'READS')
    samples = RelationshipTo('Sample', 'SAMPLE')

    @classmethod
    def from_resource(cls, assembly):
        r = Assembly(rid=assembly.id, title=assembly.name, intraspecific_name=assembly.intraspecific_name)
        r.save()
        if assembly.species_name:
            qs = Species.nodes.filter(name=assembly.species_name)
            if len(qs) == 0:
                s = Species(name=assembly.species_name)
                s.save()
            else:
                s = qs.get()

            r.species.connect(s)
        return r


class Barcodes(Resource):
    country = RelationshipTo('Country', 'COUNTRY')
    subdivision = StringProperty()
    marker = StringProperty()

    # @classmethod
    # def from_resource(cls, sample):


class Structure(Resource):
    deposit_date = DateProperty()
    method = StringProperty()

    @classmethod
    def from_resource(cls, sample):
        g = cls(rid=sample.id, title=sample.name, subdivision=sample.subdivision,
                collection_date=sample.collection_date)
        g.save()
        if sample.country:
            c = Country.nodes.get(name=sample.country)
            r.country.connect(c)
        return g


class Sample(Resource):
    country = RelationshipTo('Country', 'COUNTRY')
    subdivision = StringProperty()
    collection_date = DateProperty()

    @classmethod
    def from_resource(cls, sample):
        g = cls(rid=sample.id, title=sample.name, subdivision=sample.subdivision,
                collection_date=sample.collection_date)
        g.save()
        if sample.country:
            c = Country.nodes.get(name=sample.country)
            r.country.connect(c)
        return g


class Reads(Resource):
    sample = RelationshipTo('Sample', 'SOURCE', cardinality=One)

    @classmethod
    def from_resource(cls, reads):
        g = cls(rid=reads.id, title=reads.name)
        g.save()
        return g


class Tool(Resource):
    tool_type = StringProperty()

    @classmethod
    def from_resource(cls, tool):
        g = cls(rid=tool.id, title=tool.name, tool_type=tool.tool_type)
        g.save()
        return g


from django.db.models.signals import post_save, post_delete, post_init, m2m_changed
from django.dispatch import receiver
from bioresources.models.Resource import Collaboration
from bioresources.models.Affiliation import Affiliation

from bioresources.models.Person import Person as rPerson


@receiver(post_save, sender=Affiliation)
def affiliation_handler(sender, **kwargs):
    aff = kwargs["instance"]
    aff.author.type = rPerson.TYPE
    connect_nodes(aff.author, aff.resource)


@receiver(post_save, sender=Collaboration)
def collaboration_handler(sender, **kwargs):
    collaboration = kwargs["instance"]
    c = Collaboration.objects.prefetch_related("person", "organization", "resource").get(id=collaboration.id)
    if collaboration.person:
        c.person.type = rPerson.TYPE
        connect_nodes(c.person, c.resource, reltype=Collaboration.rev_types[c.type])
    if collaboration.organization:
        c.organization.type = rOrganization.TYPE
        connect_nodes(c.organization, c.resource, reltype=Collaboration.rev_types[c.type])




from bioresources.models.ResourceRelation import ResourceRelation


@receiver(post_save, sender=ResourceRelation)
def resource_relation_handler(sender, **kwargs):
    c = ResourceRelation.objects.prefetch_related("source", "target").get(id=kwargs["instance"].id)

    connect_nodes(c.source, c.target, reltype=c.role)


@receiver(post_delete, sender=Collaboration)
def my_handler2(sender, **kwargs):
    c = Collaboration.objects.prefetch_related("person", "resource").get(id=kwargs["instance"].id)
    c.person.type = Person.TYPE
    delete_edge(c.person, c.resource, reltype=Collaboration.rev_types[c.type])


from bioresources.models.Resource import Resource
from bioresources.models.ReadsArchive import ReadsArchive as rReadsArchive
from bioresources.models.Sample import Sample as rSample
from bioresources.models.Assembly import Assembly as rAssembly
from bioresources.models.Expression import Expression as rExpression
from bioresources.models.Barcode import Barcode as rBarcode
from bioresources.models.Structure import Structure as rStructure
from bioresources.models.Tool import Tool as rTool
from bioresources.models.Publication import Publication as rPublication
from bioresources.models.Organization import Organization as rOrganization
from bioresources.models.BioProject import BioProject as rBioProject

gclass_dict = {
    Resource.RESOURCE_TYPES.STRUCTURE: Structure,
    Resource.RESOURCE_TYPES.ASSEMBLY: Assembly,
    Resource.RESOURCE_TYPES.SAMPLE: Sample,
    Resource.RESOURCE_TYPES.EXPRESSION: Expression,
    Resource.RESOURCE_TYPES.PUBLICATION: Publication,
    Resource.RESOURCE_TYPES.TOOL: Tool,
    Resource.RESOURCE_TYPES.READS: Reads,
    Resource.RESOURCE_TYPES.BARCODE: Barcodes,

    Resource.RESOURCE_TYPES.PERSON: Person,
    Resource.RESOURCE_TYPES.ORGANIZATION: Organization,
    Resource.RESOURCE_TYPES.BIOPROJECT: BioProject

}



@receiver(post_save, sender=rBarcode)
@receiver(post_save, sender=rTool)
@receiver(post_save, sender=rReadsArchive)
@receiver(post_save, sender=rSample)
@receiver(post_save, sender=rAssembly)
@receiver(post_save, sender=rExpression)
@receiver(post_save, sender=rStructure)
@receiver(post_save, sender=rPublication)
@receiver(post_save, sender=rOrganization)
@receiver(post_save, sender=rPerson)
@receiver(post_save, sender=rBioProject)
def my_handler3(sender, **kwargs):
    r = kwargs["instance"]
    print(r)
    r = sender.objects.get(id=r.id)
    gclass = gclass_dict[sender.TYPE]
    print(gclass)
    if len(gclass.nodes.filter(rid=r.id)) == 0:
        print("xxx")
        gclass.from_resource(r)


@receiver(post_delete, sender=rBarcode)
@receiver(post_delete, sender=rTool)
@receiver(post_delete, sender=rReadsArchive)
@receiver(post_delete, sender=rSample)
@receiver(post_delete, sender=rAssembly)
@receiver(post_delete, sender=rExpression)
@receiver(post_delete, sender=rStructure)
@receiver(post_delete, sender=rPublication)
@receiver(post_delete, sender=rOrganization)
@receiver(post_delete, sender=rPerson)
@receiver(post_delete, sender=rBioProject)
def model_delete_handler(sender, **kwargs):
    r = kwargs["instance"]
    gclass = gclass_dict[sender.TYPE]
    qs = gclass.nodes.filter(rid=r.id)
    if len(qs):
        for x in qs:
            x.delete()


@receiver(m2m_changed, sender=Affiliation.organizations.through)
def aff_handler(sender, **kwargs):
    ids = kwargs.pop('pk_set')
    organization = rOrganization.objects.filter(id__in=ids)
    action = kwargs.pop('action', None)
    aff = kwargs["instance"]
    aff.author.type = rPerson.TYPE
    for org in organization:
        org.type = rOrganization.TYPE
        if action == "post_add":
            connect_nodes(aff.author, org, reltype="AFF")
        elif action in ["post_remove", "post_clear"]:
            delete_edge(aff.author, org, reltype="AFF")
