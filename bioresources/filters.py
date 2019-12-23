# import django_filters
# from django_filters.widgets import RangeWidget
from datetime import datetime

now = datetime.now()

# path('publications', views.filteredpersonlistview.as_view(), name='publications'),

def _truncate(dt):
    return dt.date()


# class PublicationFilter(django_filters.FilterSet):
#     # date_of_publication = django_filters.RangeFilter(name='date_of_publication', lookup_expr='range')
#     start_date = django_filters.DateFilter(name="date_of_publication", lookup_expr='gte')
#     end_date = django_filters.DateFilter(name="date_of_publication", lookup_expr='lte')
#
#     class Meta:
#         model = Publication
#         form = PublicationForm
#         fields = ["name", "type"]
#
#         # fields = {'name': ['contains'],
#         #           #'description': ['contains'],
#         #           date_of_publication:["c"],
#         #           "type": ["exact"]}
#
#         # filter_overrides = {
#         #
#         #     DateField: {
#         #         'filter_class': django_filters.RangeFilter,
#         #         'extra': lambda f: {
#         #             'widget': RangeWidget,
#         #         },
#         #     },
#         # }
