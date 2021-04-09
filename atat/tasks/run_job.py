# Here lies 3rd party module imports
from celery import chord, group, chain

# Here lies the Celery app import
from .atat import ATAT
from ..celery_app import app

# Here lies the model imports
from ..models import Job

# Here lies the task imports
from .hunana import hunana
from .warehousing import Warehousing
from .atat import noop

# Here lies the constant imports
from .constants import LogContexts

# Here lies other imports
from .logging import Logging


@app.task(name="Job")
def run_job(source_seqs: str, reservoir_seqs: str, jobid: str, **kwargs):
    Logging.make_log_entry(jobid, LogContexts.INFO, f'Job {jobid} is starting.', 'RUNNING')

    all_hunana_tasks = group([
        hunana.s(source_seqs, jobid, **kwargs),
        hunana.s(reservoir_seqs, jobid, **kwargs)
    ])

    atat_task = group(noop.s(), ATAT().s(jobid))

    task_chain = chain(all_hunana_tasks, atat_task)

    chord([task_chain])(Warehousing().s(jobid))

    # chain(all_hunana_tasks, ATAT2().s()).delay()


    # chord([
    #     hunana.s(Tags.SOURCE, source_seqs, **kwargs),
    #     hunana.s(Tags.RESERVOIR, reservoir_seqs, **kwargs)
    # ])(Warehousing().s(jobid)).forget()


