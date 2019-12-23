from django.conf import settings # import the settings file

def config_variables(request):

    return {'GOOGLE_ANALYTICS_CODE': settings.GOOGLE_ANALYTICS_CODE}