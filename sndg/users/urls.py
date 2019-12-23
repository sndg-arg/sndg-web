from django.urls import path

from sndg.users.views import (
    user_redirect_view,
    user_update_view,
    user_detail_view,
    task_example,
    task_res
)

from bioresources.views.user.UserResourcesView import UserResourcesView

app_name = "users"

urlpatterns = [
    path("~redirect/", view=user_redirect_view, name="redirect"),
    path("~update/", view=user_update_view, name="update"),
    path("job/", view=task_example, name="jobtest"),
    path("job/<str:jid>/", view=task_res, name="jobres"),
    path("<str:username>/", view=UserResourcesView, name="detail"),

]
