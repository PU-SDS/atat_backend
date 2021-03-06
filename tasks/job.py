# Here lies 3rd party module imports
from celery import chord

# Here lies the Celery app import
from ..celery_app import app

# Here lies the model imports
from ..models import Job

# Here lies the task imports
from .hunana import hunana
from .warehousing import Warehousing

# Here lies the constant imports
from .constants import LogContexts, Tags

# Here lies other imports
from .logging import Logging


@app.task(name="ATAT")
def atat(source_seqs: str, reservoir_seqs: str, jobid: str, **kwargs):
    Job(
        _id=jobid,
        log=[Logging.make_log_entry(context=LogContexts.INFO, msg=f'Starting job {jobid}.')],
        status='STARTED').save()

    chord([
        hunana.s(Tags.SOURCE, source_seqs, **kwargs),
        hunana.s(Tags.RESERVOIR, reservoir_seqs, **kwargs)
    ])(Warehousing().s(jobid)).forget()


