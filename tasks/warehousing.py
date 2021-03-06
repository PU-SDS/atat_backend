# Here lies the Celery app import
from celery import Task

from atat_single.celery_app import app

# Here lies the constants import
from .constants import Tags

# Here lies the model imports
from ..models import HunanaPosition, Job, Result

# Here lies other imports
from .logging import Logging

# Here lies constants imports
from .constants import LogContexts


class Warehousing(Task):
    name = "Warehousing"

    def run(self, hunana_results: list, jobid: str):
        """
            Stores the Hunana results for Source and Reservoir under the appropriate job id.

            :param hunana_results: A list containing two dicts pertaining to the tagged Hunana results for Source a
            Reservoir sequences.
            :param jobid: The current job id.

            :type hunana_results: list
            :type jobid: str
        """

        # First we get the job that we just saved using the job id. Then we update the log
        job = Job.objects.get(_id=jobid)
        job.log.append(Logging.make_log_entry(LogContexts.INFO, f'Hunana completed for job {jobid}.'))
        job.save()

        # Then we create an instance of the Results model that we will later link to the job
        result = Result()

        # The first argument to this function will always be a list containing the GROUPED Hunana results for both
        # Source and Reservoir. We loop through thatl ist, check the tag and create a list of HunanaPosition instances
        # which we will later link with the Results instance we just created above.
        for hunana_result in hunana_results:  # type: dict
            if hunana_result.get('tag') == Tags.SOURCE:
                result.source = self._get_hunana_positions(hunana_result.get('results'))
                job.log.append(
                    Logging.make_log_entry(LogContexts.INFO, f'Source sequence data stored successfully for job '
                                                             f'{jobid}'))
                job.save()
            elif hunana_result.get('tag') == Tags.RESERVOIR:
                result.reservoir = self._get_hunana_positions(hunana_result.get('results'))
                job.log.append(
                    Logging.make_log_entry(LogContexts.INFO, f'Reservoir sequence data stored successfully for '
                                                             f'job '
                                                             f'{jobid}'))
                job.save()
            else:
                return

        # We have everything we need. Link the results to the job and save. Then save the job itself
        job.results = result.save()
        job.log.append(Logging.make_log_entry(LogContexts.INFO, f'Job {jobid} successfully completed.'))
        job.save()

    @classmethod
    def _get_hunana_positions(cls, position: dict) -> list:
        position_data = [
            HunanaPosition(
                position=position.get('position'),
                entropy=position.get('entropy'),
                supports=position.get('supports'),
                variants=position.get('sequences'),
                kmer_types=position.get('kmer_types')
            )
            for position in position]

        return position_data


# Here we register the task under our  Celery app  because when ug Class based tasks auto-register does not happen
app.register_task(Warehousing())
