# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models

class ProcessStatus(models.Model):
    name = models.CharField(max_length=200, help_text=__('Process name'))
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "ProcessStatus(%s)" % self.name

    def step(self, key):
        for step in self.steps.all():
            if step.name == key:
                return step
        raise IndexError("'%s' step not found" % key)

    def __repr__(self):
        return self.__str__()


class ProcessStatusStep(models.Model):
    name = models.CharField(max_length=200)
    class_name = models.CharField(max_length=200)
    class_identifier = models.CharField(max_length=200)
    cast_int_identifier = models.BooleanField(default=True)
    created_at = models.DateTimeField(auto_now_add=True)
    completed = models.BooleanField(default=False)
    process_status = models.ForeignKey(ProcessStatus, related_name="steps", on_delete=models.CASCADE)

    def __contains__(self, key):
        return len(self.units.filter(process_identifier=key))

    def append(self, db_identifier, process_identifier=None):
        if not process_identifier:
            process_identifier = db_identifier
        ProcessStatusStepProcessUnit(process_status_step=self,
                                     db_identifier=db_identifier, process_identifier=process_identifier).save()

    def results(self):
        clazz = get_class(self.class_name)
        for unit in self.units.all():
            if unit.db_identifier:
                identifier = int(unit.db_identifier) if self.cast_int_identifier else unit.db_identifier
                yield clazz.objects.get(**{self.class_identifier: identifier})

    def __str__(self):
        return "ProcessStatusStep(%s)" % self.name

    def __repr__(self):
        return self.__str__()


class ProcessStatusStepProcessUnit(models.Model):
    db_identifier = models.CharField(max_length=30, null=True)
    process_identifier = models.CharField(max_length=50)
    created_at = models.DateTimeField(auto_now_add=True)
    process_status_step = models.ForeignKey(ProcessStatusStep, related_name="units", on_delete=models.CASCADE)

    def __str__(self):
        return "PSStepPUnit('%s','%s')" % (self.db_identifier, self.process_identifier)

    def __repr__(self):
        return self.__str__()