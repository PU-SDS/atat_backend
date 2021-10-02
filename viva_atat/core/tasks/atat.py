from typing import List

from ..models import JobDBModel, LogMessageFlags, LogMessages, JobStatus
from ...celery_app import app
from ..analysis import Analyses


@app.task(name="atat")
def atat(dima_results: List[List[dict]], job_id: str):
    """
    Runs analysis to figure out the transmission dynamics of motifs between the two datasets. This is the core
    analysis of ATAT.

    :param dima_results: This list contains two elements, the first one being the DiMA results for the source
        dataset and the second being for the reservoir (same order as Celery group).
    :param job_id: The job id we are processing so that we can update the log.

    :type dima_results: List[List[dict]]
    :type job_id: str
    """

    job_queryset = JobDBModel.objects.filter(id=job_id)
    switches = None

    try:
        job_queryset.update_log(LogMessageFlags.INFO, LogMessages.RUN_TRANSMISSIBILITY_ANALYSIS)
        switches = Analyses.at_analysis(dima_results[0], dima_results[1])
    except Exception:
        job_queryset.update_log(LogMessageFlags.ERROR, LogMessages.TRANSMISSIBILITY_ANALYSIS_ERROR)
        job_queryset.update_status(JobStatus.FAILED)
        exit(1)

    return switches


@app.task(name="noop")
def noop(hunana_results: list) -> list:
    """
    This performs no actions. Just returns the Hunana results so we can inject those results into the group
    """

    return hunana_results
