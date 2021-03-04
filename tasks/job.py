from celery.result import AsyncResult

from .hunana_reservoir import hunana_reservoir
from .hunana_source import hunana_source

from atat_single.celery_app import app


@app.task()
def run_job(source_seqs: str, reservoir_seqs: str, jobid: str, **kwargs):
    source_async = hunana_source.delay(source_seqs, **kwargs)  # type: AsyncResult
    reservoir_async = hunana_reservoir.delay(reservoir_seqs, **kwargs)  # type: AsyncResult

    source_results = source_async.get()
    reservoir_results = reservoir_async.get()

    print(reservoir_results)
