# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from model_utils import Choices
from django.utils.translation import gettext_lazy as _, ngettext as __
from django.db import models

from .Resource import Resource

class Tool(Resource):
    TYPE = Resource.RESOURCE_TYPES.TOOL

    TYPES = Choices(
        *[(i, x, _(x)) for i, x in enumerate([
            'app', 'database', 'library', 'plugin', 'program', 'webserver',
        ])]
    )

    tool_type = models.PositiveSmallIntegerField(choices=TYPES)
    url = models.URLField(null=True)

    class Meta:
        verbose_name_plural = _("Tools")