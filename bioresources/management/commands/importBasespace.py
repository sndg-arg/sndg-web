import os
import subprocess as sp

import requests
from django.core.management.base import BaseCommand
from tqdm import tqdm


class Command(BaseCommand):
    """
    https://developer.basespace.illumina.com/docs/content/documentation/rest-api/

    https://api.basespace.illumina.com/v1pre3/projects/84305427/samples?access_token=?&limit10
    https://api.basespace.illumina.com/v1pre3/samples/150888752/files?access_token=
    https://api.basespace.illumina.com/v1pre3/files/11596527869/content?access_token=
    """
    BASEPASE_URL_TEMPLATE = "https://api.basespace.illumina.com/"
    SAMPLES_PATH = "v1pre3/projects/{project_id}/samples"
    SAMPLE_FILES_PATH = "/files"
    FILE_CONTENT_PATH = "/content"
    params = {"limit": "{limit}", "offset": "{offset}", "access_token": "{access_token}"}

    help = 'Downloads sample files from basespace'

    def add_arguments(self, parser):
        parser.add_argument('--access_token', required=True)
        parser.add_argument('--project', required=True)
        parser.add_argument('--dst_dir', default="./")
        parser.add_argument('--no_color', default=False)

    def download_file(self, complete_url, target, ovewrite=False, retries=3):
        if not target.strip():
            target = "./"
        if not os.path.exists(os.path.dirname(os.path.abspath(target))):
            raise Exception("%s does not exists" % os.path.dirname(target))
        if os.path.exists(target) and not ovewrite:
            raise OvewriteFileException("%s already exists" % target)

        sp.call(' wget  --timeout=20 --tries={retries} -O {target} "{url}"'.format(
            url=complete_url, retries=retries, target=target), shell=True)

    def handle(self, *args, **options):
        assert os.path.exists(options["dst_dir"])

        options["dst_dir"] = os.path.abspath(options["dst_dir"])
        url = self.BASEPASE_URL_TEMPLATE + self.SAMPLES_PATH + "?" + "&".join(
            [k + "=" + v for k, v in self.params.items()])
        r = requests.get(
            url.format(access_token=options["access_token"],
                       limit=1000, offset=0, project_id=options["project"]))
        samples = []
        for item in r.json()["Response"]["Items"]:
            samples.append(item['Href'])
        for sample in tqdm(samples):
            url = self.BASEPASE_URL_TEMPLATE + sample + self.SAMPLE_FILES_PATH + "?" + "&".join(
                [k + "=" + v for k, v in self.params.items()])
            r = requests.get(url.format(access_token=options["access_token"],
                                        limit=1000, offset=0))
            for item in r.json()["Response"]["Items"]:
                url = self.BASEPASE_URL_TEMPLATE + item["HrefContent"] + "?access_token=" + options["access_token"]
                self.download_file(url, options["dst_dir"] + "/" + item["Name"])
