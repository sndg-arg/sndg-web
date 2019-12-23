from django.urls import path

from . import views

app_name = 'vardb'
urlpatterns = [
    path('variant/<int:pk>', views.variant, name='variant_view'),
    path('variant_collection/<int:pk>', views.variant_collection, name='variant_collection_view'),
    path('collection_set/<int:pk>', views.collection_set, name='collection_set_view'),



]