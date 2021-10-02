from celery import chord, group, chain

from .atat import atat, noop
from .warehousing import warehousing
from ..models import JobDBModel, LogMessageFlags, JobStatus, LogMessages
from ...celery_app import app

from .dima import dima


@app.task(name="job")
def run_job(source_seqs: str, reservoir_seqs: str, job_id: str):
    # We get the queryset so we can update the logs easier
    job_queryset = JobDBModel.objects.filter(id=job_id)

    # We then retrieve the job from the database
    job = job_queryset.get()

    # Update the log to indicate that we have started
    job_queryset.update_log(LogMessageFlags.INFO, LogMessages.JOB_RUNNING)
    job_queryset.update_status(JobStatus.RUNNING)

    all_dima_tasks = group(
        [
            dima.s(source_seqs, job_id, job.parameters.kmer_length, job.parameters.header_format),
            dima.s(reservoir_seqs, job_id, job.parameters.kmer_length, job.parameters.header_format),
        ]
    )

    atat_task = group(noop.s(), atat.s(job_id))

    task_chain = chain(all_dima_tasks, atat_task)

    chord([task_chain])(warehousing.s(job_id))
