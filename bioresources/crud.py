from crudbuilder.abstract import BaseCrudBuilder

from .models import Publication


class PublicationCrud(BaseCrudBuilder):
    model = Publication
    search_fields = ['name']
    tables2_fields = ('name', 'description')
    tables2_css_class = "table table-bordered table-condensed"
    tables2_pagination = 20  # default is 10
    #modelform_excludes = ['created_by', 'updated_by']
    login_required=True
    permission_required=True

    # @classmethod
    # def custom_queryset(cls, request, **kwargs):
    #     """Define your own custom queryset for list view"""
    #     qset = cls.model.objects.filter(created_by=request.user)
    #     return qset
    #
    # @classmethod
    # def custom_context(cls, request, context, **kwargs):
    #     """Define your own custom context for list view"""
    #     context['custom_data'] = "Some custom data"
    #     return context

    permissions = {
        'list': 'example.person_list',
        'create': 'example.person_create'
    }
    # createupdate_forms = {
    #     'create': PersonCreateForm,
    #     'update': PersonUpdateForm
    # }