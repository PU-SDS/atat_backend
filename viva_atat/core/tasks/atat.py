from typing import List

from ..models import JobDBModel, LogMessageFlags, LogMessages, JobStatus, Transmission, Results
from ...celery_app import app
from ..analysis import Analyses


@app.task(name="atat-standalone")
def atat_standalone(dima_results: List[List[dict]], job_id: str):
    """
    Runs analysis to figure out the transmission dynamics of motifs between the two datasets. This is the core
    analysis of ATAT. This is the task that runs when ATAT is running in standalone mode.

    :param dima_results: This list contains two elements, the first one being the DiMA results for the host
        dataset and the second being for the reservoir (same order as Celery group).
    :param job_id: The job id we are processing so that we can update the log.

    :type dima_results: List[List[dict]]
    :type job_id: str

    :returns: A list of motif transmissions between the two datasets.
    """

    job_queryset = JobDBModel.objects.filter(id=job_id)
    transmissions = None

    try:
        job_queryset.update_log(LogMessageFlags.INFO, LogMessages.RUN_TRANSMISSIBILITY_ANALYSIS)
        transmissions = Analyses.at_analysis(dima_results[0], dima_results[1])
    except Exception:
        job_queryset.update_log(LogMessageFlags.ERROR, LogMessages.TRANSMISSIBILITY_ANALYSIS_ERROR)
        job_queryset.update_status(JobStatus.FAILED)
        exit(1)

    return transmissions


@app.task(name="atat-viva")
def atat_viva(job_id: str):
    """
    Runs analysis to figure out the transmission dynamics of motifs between the two datasets. This is the core
    analysis of ATAT. This is the task that runs when ATAT is running as part of ViVA, where the DiMA results are
    provided by DIM (Data Integration Module) thereby eliminating the need to run DiMA ourselves.

    The job is registered already at the endpoint, here we just pull the host and reservoir dataset using the job id
    and carry on.

    :param job_id: The job id we are processing so that we can update the log.
    :type job_id: str
    """

    job_queryset = JobDBModel.objects.filter(id=job_id)
    job = job_queryset.get()  # type: JobDBModel
    transmissions = None

    job_queryset.update_status(JobStatus.RUNNING)

    host_dima_positions = [position.to_mongo().to_dict() for position in job.results.host]
    reservoir_dima_positions = [position.to_mongo().to_dict() for position in job.results.reservoir]

    try:
        job_queryset.update_log(LogMessageFlags.INFO, LogMessages.RUN_TRANSMISSIBILITY_ANALYSIS)
        transmissions = Analyses.at_analysis(host_dima_positions, reservoir_dima_positions)
        job_queryset.update_log(LogMessageFlags.INFO, LogMessages.TRANSMISSIBILITY_ANALYSIS_COMPLETE)
    except Exception:
        job_queryset.update_log(LogMessageFlags.ERROR, LogMessages.TRANSMISSIBILITY_ANALYSIS_ERROR)
        job_queryset.update_status(JobStatus.FAILED)
        exit(1)

    job_queryset.update_log(LogMessageFlags.INFO, LogMessages.SAVING_TO_DB)

    results = Results.objects.get(id=job.results.id)
    results.switches = [Transmission(**transmission.dict()) for transmission in transmissions]
    results.save()

    job_queryset.update_log(LogMessageFlags.INFO, LogMessages.ADDED_MOTIF_SWITCH_RESULTS)
    job_queryset.update_log(LogMessageFlags.SUCCESS, LogMessages.JOB_COMPLETED)
    job_queryset.update_status(JobStatus.COMPLETED)


@app.task(name="noop")
def noop(dima_results: list) -> list:
    """
    This performs no actions. Just returns the DiMA results so we can inject those results into the group
    """

    return dima_results
