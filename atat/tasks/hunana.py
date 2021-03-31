# Here lies the Python module imports
from io import StringIO

# Here lies the Celery app import
from ..celery_app import app

# Here lies 3rd party module imports
from hunana import Hunana


@app.task(name="Hunana")
def hunana(seqs: str, **kwargs) -> list:
    """
        Runs Hunana. Any kwargs will be directly passed to Hunana.

        :param seqs: The sequence to run through Hunana
        :type seqs: str

        :returns: A dictionary containing a tag and Hunana results.
    """

    results = Hunana(StringIO(seqs), **kwargs).run()

    return results
