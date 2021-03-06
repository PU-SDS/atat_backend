# Here lies the Python module imports
from io import StringIO

# Here lies the Celery app import
from atat_single.celery_app import app

# Here lies 3rd party module imports
from hunana import Hunana


@app.task(name="Hunana")
def hunana(tag: str, seqs: str, **kwargs) -> dict:
    """
        Runs Hunana.
        Any kwargs will be directly passed to Hunana.

        :param tag: The tag for the Hunana job (source/reservoir)
        :param seqs: The sequence to run through Hunana

        :type tag: str
        :type seqs: str

        :returns: A dictionary containing a tag and Hunana results.
    """

    results = Hunana(StringIO(seqs), **kwargs).run()

    return {'tag': tag, 'results': results}
