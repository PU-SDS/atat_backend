# Here lies the Celery app import
from celery import Task

from ..celery_app import app

# Here lies the model imports
from ..models import HunanaPosition, Job, Result, Switch

# Here lies other imports
from .logging import Logging

# Here lies constants imports
from .constants import LogContexts, JOB_ID_GLOBAL

# Only for typing
from .atat import ATAT
from typing import List, Tuple, Union


class Warehousing(Task):
    name = "Warehousing"
    queue = "Warehousing"

    def run(self, results: list):
        """
            Stores the Hunana results for Source and Reservoir under the appropriate job id.

            :param hunana_results: A list containing two dicts pertaining to the Hunana results for Source a
            Reservoir sequences.
            :param jobid: The current job id.

            :type hunana_results: list
            :type jobid: str
        """

        Logging.make_log_entry(LogContexts.INFO, 'Storing k-mer data in database.')

        hunana_results, atat_results = results[0]  # type: Union[list, List[ATAT.switch]]
        source_hunana_results, reservoir_hunana_results = hunana_results  # type: dict

        # First we get the job that we just saved using the job id. Then we update the log
        job = Job.objects.get(_id=JOB_ID_GLOBAL)

        # Then we create an instance of the Results model that we will later link to the job
        result = Result()

        result.source = self._get_hunana_positions(source_hunana_results)
        result.reservoir = self._get_hunana_positions(reservoir_hunana_results)
        job.save()

        Logging.make_log_entry(LogContexts.INFO, 'K-mer data storing completed.')
        Logging.make_log_entry(LogContexts.INFO, 'Storing transmissibility analysis data in database.')

        # Now we need to save the ATAT (motif switching) data
        switches = [Switch(
            position=motif_switch.get('position'),
            sequence=motif_switch.get('sequence'),
            fromx=motif_switch.get('from'),
            to=motif_switch.get('to')
        ) for motif_switch in atat_results]

        result.switches = switches

        Logging.make_log_entry(LogContexts.INFO, 'Transmissibility analysis data stored successfully.')

        # We have everything we need. Link the results to the job and save. Then save the job itself
        job.results = result.save()
        Logging.make_log_entry(LogContexts.INFO, 'Job completed successfully.')
        job.status = 'FINISHED'
        job.save()

    @classmethod
    def _get_hunana_positions(cls, position: dict) -> list:
        position_data = [
            HunanaPosition(
                position=position.get('position'),
                entropy=position.get('entropy'),
                supports=position.get('supports'),
                variants=position.get('variants')
            )
            for position in position]

        return position_data


# Here we register the task under our  Celery app  because when ug Class based tasks auto-register does not happen
app.register_task(Warehousing())
