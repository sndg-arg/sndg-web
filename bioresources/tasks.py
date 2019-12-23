# -*- coding: utf-8 -*-
from __future__ import absolute_import, unicode_literals
from django.utils.translation import gettext as _

import os
from celery import shared_task

from django.core.mail import send_mail
from django.conf import settings
from sentry_sdk import capture_message

@shared_task(soft_time_limit=60*60, time_limit=60*60+20)
def execute_job(jobId):
    from .models.Job import Job
    # if not os.path.exists("/tmp/sndg"): #TODO hacerlo configurable
    #     os.makedirs("/tmp/sndg")
    job = Job.objects.get(id=jobId)
    job.running()
    job.save()
    try:
        job.execute()
    except Exception as ex:
        job.error()
    with open(job.dev_error) as h:
        error = h.read()
        if error:
            capture_message(error)

    job.save(force_update=True)
    if job.user:
        send_mail(
            _('Job results for ID %(jid)s') % {"jid":job.id},
            _('Click <a href="%(link)s">here</a> to see results for job %(jid)s') % {
                "jid":job.id,"link":job.get_absolute_url()} ,
            settings.EMAIL_HOST_USER,
            [job.user.email],
            fail_silently=False,
        )
