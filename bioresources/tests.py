# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TestCase
from .models.Organization import Organization
from django.test import Client
csrf_client = Client(enforce_csrf_checks=True)

# from bioresources.io.adapters import NCBIGDSAdapter


class AdaptersTestCase(TestCase):
    def setUp(self):
        pass

    # def test_basic_data(self):
    #     Organization.objects.get(name="NCBI")
    #     Organization.objects.get(name="Bold")
    #     # adapter = NCBIGDSAdapter()
    #     # summaryData = adapter.fetch(GSE107376)
    #     # adapter.save(summaryData)
    #     #
    #     # self.assertEqual(lion.speak(), 'The lion says "roar"')
    #     # self.assertEqual(cat.speak(), 'The cat says "meow"')

    def test_basic_urls(self):
        for url2Test in ["/"]:
            csrf_client.get(url2Test)

