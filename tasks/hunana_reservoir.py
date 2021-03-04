from io import StringIO

from atat_single.celery_app import app
from hunana import Hunana


@app.task()
def hunana_reservoir(seqs: str, **kwargs):
    return Hunana(StringIO(seqs), **kwargs).run()
