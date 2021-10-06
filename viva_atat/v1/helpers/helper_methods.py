import json

from mongoengine import DoesNotExist

from .exceptions import JobExists
from ...core.models import JobDBModel, Parameters, Results, DimaPosition, LogMessageFlags, LogMessages
from ..models import CreateStandaloneJobRequest, CreateVivaJobRequest, CreateJobParameters
from ...core.tasks import atat_viva, run_job


class HelperMethods(object):
    @staticmethod
    def _create_db_parameter_obj(parameters: CreateJobParameters) -> Parameters:
        """
        Private method that creates a parameter object, saves it in the database and returns a Parameters object.

        :param parameters: The parameters to be used for the job.
        :type parameters: CreateJobParameters

        :returns: A DB parameters object that can be linked with a job in the database.
        """

        return Parameters(kmer_length=parameters.kmer_length, header_format=parameters.header_format).save()

    @classmethod
    def _check_job_exists(cls, job_id: str) -> bool:
        """
        Given a job id check if it exits already in the database.

        :param job_id: The job id provided by the user.
        :type job_id: str

        :returns: A boolean value indicating whether the job exists or not.
        """

        status = True

        try:
            cls.get_job(job_id)
        except DoesNotExist:
            status = False

        return status

    @classmethod
    def create_standalone_job(cls, payload: CreateStandaloneJobRequest) -> str:
        """
        Given a valid set of parameters originating from a valid standalone job create request, registers it as a new
        job on the database and returns the auto-generated job id.

        :param payload: A payload originating from a valid request to the standalone job creation endpoint.
        :type payload: CreateStandaloneJobRequest

        :returns: The auto-generated job id.
        """

        parameters = cls._create_db_parameter_obj(payload.parameters)
        job = JobDBModel(parameters=parameters).save()

        run_job.delay(payload.host_sequences, payload.reservoir_sequences, job.id)

        return job.id

    @classmethod
    def create_viva_job(cls, payload: CreateVivaJobRequest, testing: bool = False) -> str:
        """
        Given a valid set of parameters originating from a valid job create request, registers it as a new job on the
        database and returns the auto-generated job id.

        :param payload: The payload originating from a valid request to the viva job creation endpoint.
        :param testing: Whether we are running in testing mode.

        :type payload: CreateVivaJobRequest
        :type testing: bool

        :returns: The Celery task id.
        """

        # First we find out of the job exists, if so, we raise a custom error that we can handle upstream
        status = cls._check_job_exists(payload.id)
        if status:
            raise JobExists

        # Register the job parameters in the database
        parameters = cls._create_db_parameter_obj(payload.parameters)

        # Register the job and get the queryset for updating logs
        job = JobDBModel(id=payload.id, parameters=parameters).save()
        job_queryset = JobDBModel.objects.filter(id=job.id)

        # Update the log
        job_queryset.update_log(LogMessageFlags.INFO, LogMessages.JOB_CREATED)

        # Create save the DiMA results in the database
        host_positions = [DimaPosition(**json.loads(position.json())) for position in payload.host_dima_positions]
        reservoir_positions = [
            DimaPosition(**json.loads(position.json())) for position in payload.reservoir_dima_positions
        ]
        results = Results(host=host_positions, reservoir=reservoir_positions).save()

        # Save the dima results
        job.results = results
        job.save()

        # Update the log to say that we stored the dima results
        job_queryset.update_log(LogMessageFlags.INFO, LogMessages.ADDED_MOTIF_RESULTS)

        # Then we throw the ATAT task to the queue
        if testing:
            task_id = atat_viva(job.id)
        else:
            task_id = atat_viva.delay(job.id)

        return task_id

    @staticmethod
    def get_job(job_id: str) -> JobDBModel:
        """
        Given a valid job id, returns an instance of the job model pertaining to the job id.

        :param job_id: The ID of the job to retrieve.
        :type job_id: str

        :returns: An instance of a job model pertaining to the job id.
        """

        return JobDBModel.objects.get(id=job_id)

    @staticmethod
    def get_results_queryset(results_id: str):
        """
        Given a valid job id, returns the queryset of the job model pertaining to the job id.

        :param results_id: The ID of the results.
        :type results_id: str

        :returns: The queryset of a results model pertaining to the results id.
        """

        return Results.objects.filter(id=results_id)
