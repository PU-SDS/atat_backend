from ..models import (
    JobDBModel,
    LogMessageFlags,
    LogMessages,
    ResultsDBModel,
    DimaPositionDBModel,
    TransmissionDBModel,
    JobStatus,
)
from ...celery_app import app


def _generate_pos_objs(positions, job_id: str, ds: int):
    pos_objs = list()

    for position in positions:
        pos_dict = dict()

        for k, v in position.items():
            # noinspection PyProtectedMember
            if k in DimaPositionDBModel._fields.keys():
                if k == 'position':
                    v = f'{v}|{job_id}|{ds}'

                pos_dict.update({k: v})

        pos_objs.append(DimaPositionDBModel(**pos_dict).save())

    return pos_objs


@app.task(name="warehousing")
def warehousing(results: list, job_id: str):
    """
    Stores all the results generated in the database under the correct job id.

    :param results: All the results generated upstream are contained within this variable.
                |_____ dima
                |        |______ dataset one
                |        |______ dataset two
                |_____ atat
    :param job_id: The job id we are processing currently.

    :type results: list
    :type job_id: str
    """

    # We get the queryset so we can update the logs easier
    job_queryset = JobDBModel.objects.filter(id=job_id)

    # We then retrieve the job from the database
    job = job_queryset.get()

    # Update log
    job_queryset.update_log(LogMessageFlags.INFO, LogMessages.SAVING_TO_DB)

    # Separate the data from the aggregate
    dima_results, atat_results = results[0]
    dataset_one_dima_results, dataset_two_dima_results = dima_results

    # Then we create an instance of the Results model that we will later link to the job
    result = ResultsDBModel()

    # Add the kmer positions to the database
    result.dataset_one = list(_generate_pos_objs(dataset_one_dima_results, job_id, 1))
    result.dataset_two = list(_generate_pos_objs(dataset_two_dima_results, job_id, 2))

    job.save()

    job_queryset.update_log(LogMessageFlags.INFO, LogMessages.ADDED_MOTIF_RESULTS)

    # Now we need to save the ATAT (motif switching) data
    switches = [TransmissionDBModel(**motif_switch).save() for motif_switch in atat_results]

    result.switches = switches

    job_queryset.update_log(LogMessageFlags.INFO, LogMessages.ADDED_MOTIF_SWITCH_RESULTS)

    # We have everything we need. Link the results to the job and save. Then save the job itself
    job.results = result.save()

    job_queryset.update_log(LogMessageFlags.SUCCESS, LogMessages.JOB_COMPLETED)
    job_queryset.update_status(JobStatus.COMPLETED)

    job.save()
