# Here lies the Celery app import
from celery import Task

from atat_single.celery_app import app

# Here lies the model imports
from ..models import HunanaPosition, Job, Result, Switch

# Here lies other imports
from .logging import Logging

# Here lies constants imports
from .constants import LogContexts

# Only for typing
from .atat import ATAT
from typing import List, Tuple, Union


class Warehousing(Task):
    name = "Warehousing"
    queue = "Warehousing"

    def run(self, results: list, jobid: str):
        """
            Stores the Hunana results for Source and Reservoir under the appropriate job id.

            :param hunana_results: A list containing two dicts pertaining to the Hunana results for Source a
            Reservoir sequences.
            :param jobid: The current job id.

            :type hunana_results: list
            :type jobid: str
        """

        hunana_results, atat_results = results[0]  # type: Union[list, List[ATAT.switch]]
        source_hunana_results, rervoir_hunana_results = hunana_results  # type: dict

        # First we get the job that we just saved using the job id. Then we update the log
        job = Job.objects.get(_id=jobid)
        job.log.append(Logging.make_log_entry(LogContexts.INFO, f'Hunana completed for job {jobid}.'))
        job.status = 'RUNNING'
        job.save()

        # Then we create an instance of the Results model that we will later link to the job
        result = Result()

        result.source = self._get_hunana_positions(source_hunana_results)
        result.reservoir = self._get_hunana_positions(rervoir_hunana_results)

        job.log.append(
            Logging.make_log_entry(LogContexts.INFO, f'Reservoir sequence data stored successfully for '
                                                     f'job '
                                                     f'{jobid}'))
        job.save()

        # Now we need to save the ATAT (motif switching) data
        print(atat_results)
        switches = [Switch(
            position=motif_switch.get('position'),
            sequence=motif_switch.get('sequence'),
            fromx=motif_switch.get('from'),
            to=motif_switch.get('to')
        ) for motif_switch in atat_results]

        result.switches = switches

        # We have everything we need. Link the results to the job and save. Then save the job itself
        job.results = result.save()
        job.log.append(Logging.make_log_entry(LogContexts.INFO, f'Job {jobid} successfully completed.'))
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
