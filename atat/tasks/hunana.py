# Here lies the Python module imports
from io import StringIO

# Here lies the Celery app import
from .logging import Logging
from .run_job import JOB_ID_GLOBAL
from ..celery_app import app

# Here lies 3rd party module imports
from hunana import Hunana

from ..models import Job


@app.task(name="Hunana")
def hunana(seqs: str, **kwargs) -> list:
    """
        Runs Hunana. Any kwargs will be directly passed to Hunana.

        :param seqs: The sequence to run through Hunana
        :type seqs: str

        :returns: A dictionary containing a tag and Hunana results.
    """

    Job(_id=JOB_ID_GLOBAL).log.append(
        Logging.make_log_entry('INFO', 'Generating k-mers.')
    )

    results = Hunana(StringIO(seqs), **kwargs).run()

    Job(_id=JOB_ID_GLOBAL).log.append(
        Logging.make_log_entry('INFO', 'Generating k-mers completed.')
    )

    return results
