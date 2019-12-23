# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.shortcuts import redirect, reverse, render

from django.core.exceptions import ImproperlyConfigured
from django.core.files.storage import FileSystemStorage
from django.http import HttpResponse, HttpResponseRedirect
from django import forms
from django.db import transaction

from bioresources.models.Resource import Resource
from bioresources.models.jobs.LoadGenomeJob import LoadGenomeFromFileJob
from bioseq.models.Biodatabase import Biodatabase

def submission(request):
    return render(request, "forms/submission.html")


from django.views import generic


class UploadView(generic.TemplateView):
    template_name = 'forms/resource_upload.html'


from django.conf import settings

from resumable.fields import ResumableFileField


class ResumableForm(forms.Form):
    file = ResumableFileField(
        # allowed_mimes=("audio/ogg",),
        upload_url=lambda: reverse('upload_api'),
        chunks_dir=getattr(settings, 'FILE_UPLOAD_TEMP_DIR')
    )


from resumable.files import ResumableFile
from bioresources.tasks import execute_job

def upload_view(request):
    form = ResumableForm()
    return render(request, 'forms/resource_upload.html', {'form': form})


class ResumableUploadView(generic.TemplateView):
    template_name = 'forms/resource_upload.html'

    def get(self, request, *args, **kwargs):


        if "resumableChunkNumber" in request.GET:
            r = ResumableFile(self.storage, self.request.GET)
            if not (r.chunk_exists or r.is_complete):
                return HttpResponse('chunk not found', status=404)
            return HttpResponse('chunk already exists')
        else:
            resource = Resource.objects.get(id=kwargs["resource_id"])
            bdb = Biodatabase.objects.filter(name=resource.name)
            loaded = bool(bdb.count())
            if loaded:
                return HttpResponseRedirect(resource.get_absolute_url())

            form = ResumableForm()
            return render(request, self.template_name, {"resource": resource,
                                                        'form': form, "resource_id": kwargs["resource_id"]})

    def post(self, *args, **kwargs):
        """Saves chunks then checks if the file is complete.
        """
        chunk = self.request.FILES.get('file')
        r = ResumableFile(self.storage, self.request.POST)
        if r.chunk_exists:
            return HttpResponse('chunk already exists')
        r.process_chunk(chunk)
        if r.is_complete:
            resource = Resource.objects.get(id=kwargs["resource_id"])
            self.process_file(r.filename, r,resource)
            r.delete_chunks()

            return HttpResponseRedirect(resource.get_absolute_url())
        return HttpResponse()

    def process_file(self, filename, resumablefile,resource):
        """ Process the complete file. """
        self.storage.save(filename, resumablefile)
        with transaction.atomic():
            job = LoadGenomeFromFileJob(assembly=resource,
                                        filename=settings.FILE_UPLOAD_TEMP_DIR + "/" + filename)
            job.save()
            job.init()
            job.queue()
            job.save()

        execute_job.apply_async(args=(job.id,),countdown=10)

    @property
    def chunks_dir(self):
        chunks_dir = getattr(settings, 'FILE_UPLOAD_TEMP_DIR', None)
        if not chunks_dir:
            raise ImproperlyConfigured(
                'You must set settings.FILE_UPLOAD_TEMP_DIR')
        return chunks_dir

    @property
    def storage(self):
        return FileSystemStorage(location=self.chunks_dir)

# from django_fine_uploader.views import FineUploaderView
# class NotConcurrentUploaderView(FineUploaderView):
#     """Example of a chunked, but NOT concurrent upload.
#     Disabling concurrent uploads per view.
#     Remember, you can turn off concurrent uploads on your settings, with:
#     FU_CONCURRENT_UPLOADS = False
#     """
#     @property
#     def concurrent(self):
#         return False
#
#     def form_valid(self, form):
#         self.process_upload(form)
#         return self.make_response({'success': True})
#

# path('upload/', view=views.ExampleView.as_view(), name='home'),
# path('upload_api/', view=views.NotConcurrentUploaderView.as_view(), name='upload'),
# path('upload_api/', view=ResumableUploadView.as_view(),  name='upload_api'),
