# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.conf import settings
from django.shortcuts import render, redirect, reverse
from django.db import transaction
from datetime import datetime, timezone

from bioresources.models.Job import Job,CmdJob
from bioresources.models.BlastDB import BlastDB
from bioresources.tasks import execute_job

import Bio.SeqIO as bpio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq



def blast(request):

    if request.method == 'POST':
        """
        


'evalue_select' (139718239766000) = {str} '10'
'low_complexity_select' (139718239744696) = {str} 'Yes'
'results_select' (139718239766064) = {str} '10'
        """

        cmd = request.POST["blastType"]
        job = CmdJob(command=cmd,result_type="blast")
        if not request.user.is_anonymous:
            job.user=request.user
        job.save()
        job.init()
        job.save()
        query = settings.JOBSDIR + "/" + str(job.id) + "/query.fasta"
        with open(query,"w") as h:
            seq = "".join(request.POST["sequence_text"].split())
            r = SeqRecord(id="query",name="",description="",seq=Seq(seq))
            bpio.write(r,h,"fasta")
        db = settings.BLASTDBSDIR + "/" +  request.POST["database_select"] + ".fasta"
        kwargs = {
            "db":db,
            "query":query,
            #"matrix":request.POST["matrix_select"],
            "gapopen":request.POST["gap_open_input"],
            "gapextend":request.POST["gap_extend_input"],
            "evalue":request.POST["evalue_select"],
            "soft_masking":request.POST["low_complexity_select"],
            "max_target_seqs":request.POST["results_select"],
            "outfmt":"5"

        }
        for k,v in kwargs.items():
            cmd += " -" + k + " " + v
        job.command = cmd
        with transaction.atomic():
            job.queue()
            job.save()

        execute_job.apply_async(args=(job.id,),countdown=10)


        return redirect(reverse("bioresources:job_view", kwargs={"jid": job.id}))

    return render(request, 'tools/blast.html',
                  {"ndatabases": BlastDB.objects.filter(dbtype="nucl"),
                   "pdatabases": BlastDB.objects.filter(dbtype="prot")
                   })


def job_view(request, jid):
    job = Job.objects.get(id=jid)
    if job.status == Job.STATUS.FINISHED:
        with open(job.result) as h:
            result = h.read() #.replace("\n"," ")


        return render(request,   'tools/%s_result.html' % job.result_type, {"job":job,"result":result})
    else:
        elapsed_time = int((datetime.now(timezone.utc) - job.start).total_seconds())
        return render(request, 'tools/generic_result.html', {"job":job,"status":str(Job.STATUS[job.status]),
                                                             "elapsed_time": elapsed_time})
