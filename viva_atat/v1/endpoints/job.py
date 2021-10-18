from fastapi import APIRouter, status, HTTPException
from mongoengine import DoesNotExist

from ..helpers.exceptions import JobExists
from ..models import CreateStandaloneJobRequest, JobLogs, JobLogEntry, Parameters, CreateVivaJobRequest
from ..helpers import HelperMethods
from ...core.models import JobStatus

router = APIRouter(prefix='/job', tags=['job'])


@router.post(
    '/create/standalone', status_code=status.HTTP_201_CREATED, response_description="Returns the auto-generated job id."
)
def create_standalone_job(payload: CreateStandaloneJobRequest) -> str:
    """
    Submit a new standalone analysis job.

    A successful request will return a HTTP 201 status code, while any invalid request (missing fields) will return
    an HTTP 500 status code with the body indicating the missing fields.
    """

    job_id = HelperMethods.create_standalone_job(payload)

    return job_id


@router.post('/create/viva', status_code=status.HTTP_201_CREATED, response_description="Returns the Celery task id")
def create_viva_job(payload: CreateVivaJobRequest) -> str:
    """
    Submit a new ViVA analysis job.

    A successful request will return a HTTP 201 status code, and the body will contain the Celery task ID while any
    invalid request (missing fields) will return an HTTP 500 status code with the body indicating the missing fields.

    If the job id already exists in the ATAT database, an HTTP 409 status will be returned.
    If some unknown error occurs, an HTTP 500 error will be raised with the body indicating the reason.

    **If the position objects contain extra keys that are not required for ATAT, they will be ignored by ATAT.**
    Therefore it is safe to pass the DiMA positions as is without removing the extra fields.
    """

    task_id = None

    try:
        task_id = HelperMethods.create_viva_job(payload)
    except JobExists:
        raise HTTPException(status_code=status.HTTP_409_CONFLICT, detail="Job id already exists.")
    except Exception:
        HelperMethods.delete_job(payload.id)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Job creation failed. You may try again with the same job id.",
        )

    return str(task_id)


@router.get(
    '/status/{job_id}',
    status_code=status.HTTP_200_OK,
    response_description="Returns the current status of the job",
    response_model=JobStatus,
)
def get_job_status(job_id: str) -> JobStatus:
    """
    Get the status of a job.

    If the provided job id is not found an HTTP 404 status is returned. A successful request will return an HTTP 200.
    """

    try:
        job = HelperMethods.get_job(job_id)
    except DoesNotExist:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Job not found, consider creating one.")

    return job.status


@router.get(
    '/log/{job_id}',
    status_code=status.HTTP_200_OK,
    response_description="Returns the current status of the job",
    response_model=JobLogs,
)
def get_job_log(job_id: str) -> JobLogs:
    """
    Get the logs of a job.

    If the provided job id is not found an HTTP 404 status is returned. A successful request will return an HTTP 200.
    """

    try:
        job = HelperMethods.get_job(job_id)
    except DoesNotExist:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Job not found, consider creating one.")

    entries = [JobLogEntry(**entry.to_mongo().to_dict()) for entry in job.log]

    return JobLogs(logs=entries)


@router.get(
    '/parameters/{job_id}',
    status_code=status.HTTP_200_OK,
    response_description="Returns the parameters that were used to run the job",
    response_model=Parameters,
)
def get_job_parameters(job_id: str) -> Parameters:
    """
    Get the parameters that were used to run the job.

    If the provided job id is not found an HTTP 404 status is returned. A successful request will return an HTTP 200.
    """

    try:
        job = HelperMethods.get_job(job_id)
    except DoesNotExist:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Job not found, consider creating one.")

    return Parameters(**job.parameters.to_mongo().to_dict())
