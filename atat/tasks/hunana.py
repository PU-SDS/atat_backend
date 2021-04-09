# Here lies the Python module imports
from io import StringIO

# Here lies the Celery app import
from .logging import Logging
from .constants import JOB_ID_GLOBAL, LogContexts
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

    Logging.make_log_entry(LogContexts.INFO, 'Generating k-mers using HUNANA algorithm.')

    results = Hunana(StringIO(seqs), **kwargs).run()

    Logging.make_log_entry(LogContexts.INFO, 'K-mer generation completed successfully.')

    return results
