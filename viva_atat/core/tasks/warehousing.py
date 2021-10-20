from ..models import (
    JobDBModel,
    LogMessageFlags,
    LogMessages,
    Results,
    DimaPosition,
    Transmission,
    JobStatus,
)
from ...celery_app import app


def _generate_pos_objs(positions):
    for position in positions:
        # noinspection PyProtectedMember
        filtered = {k: v for k, v in position.items() if k in DimaPosition._fields.keys()}
        yield DimaPosition(**filtered)


@app.task(name="warehousing")
def warehousing(results: list, job_id: str):
    """
    Stores all the results generated in the database under the correct job id.

    :param results: All the results generated upstream are contained within this variable.
                |_____ dima
                |        |______ host
                |        |______ reservoir
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
    host_dima_results, reservoir_dima_results = dima_results

    # Then we create an instance of the Results model that we will later link to the job
    result = Results()

    # Add the kmer positions to the database
    result.host = list(_generate_pos_objs(host_dima_results))
    result.reservoir = list(_generate_pos_objs(reservoir_dima_results))

    job.save()

    job_queryset.update_log(LogMessageFlags.INFO, LogMessages.ADDED_MOTIF_RESULTS)

    # Now we need to save the ATAT (motif switching) data
    switches = [Transmission(**motif_switch) for motif_switch in atat_results]

    result.switches = switches

    job_queryset.update_log(LogMessageFlags.INFO, LogMessages.ADDED_MOTIF_SWITCH_RESULTS)

    # We have everything we need. Link the results to the job and save. Then save the job itself
    job.results = result.save()

    job_queryset.update_log(LogMessageFlags.SUCCESS, LogMessages.JOB_COMPLETED)
    job_queryset.update_status(JobStatus.COMPLETED)

    job.save()
