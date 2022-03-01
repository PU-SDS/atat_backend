from mongoengine import DoesNotExist, LazyReferenceField

from .exceptions import JobExists
from ...core.models import JobDBModel, ParametersDBModel, ResultsDBModel, DimaPositionDBModel, LogMessageFlags, LogMessages
from ..models import CreateStandaloneJobRequest, CreateVivaJobRequest, Parameters
from ...core.tasks import atat_viva, run_job


class HelperMethods(object):
    @staticmethod
    def _create_db_parameter_obj(parameters: Parameters) -> ParametersDBModel:
        """
        Private method that creates a parameter object, saves it in the database and returns a Parameters object.

        :param parameters: The parameters to be used for the job.
        :type parameters: Parameters

        :returns: A DB parameters object that can be linked with a job in the database.
        """

        return ParametersDBModel(kmer_length=parameters.kmer_length, header_format=parameters.header_format,
                                 protein_name=parameters.protein_name).save()

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

        run_job.delay(payload.dataset_one, payload.dataset_two, job.id)

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
        dataset_one_positions = [DimaPositionDBModel(**position.dict()) for position in payload.dataset_one_positions]
        dataset_two_positions = [
            DimaPositionDBModel(**position.dict()) for position in payload.dataset_two_positions
        ]
        results = ResultsDBModel(dataset_one=dataset_one_positions, dataset_two=dataset_two_positions).save()

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

        return ResultsDBModel.objects.filter(id=results_id)

    @staticmethod
    def delete_job(job_id: str):
        """
        Deletes a job given the job id.

        :param job_id: The ID of the job to delete.
        :type job_id: str
        """

        HelperMethods.get_job(job_id).delete()

    @staticmethod
    def get_id_from_pk(pk: str) -> int:
        """
        Since we are using a whole new document for each dima k-mer position, it gets expensive when we need to get
        a single position by position number (because have to dereference each document).
        A solution to this is to use LazyReference field which gives us access to the pk without dereferencing.
        So we store the position in the pk like so "1|jobid". We then split it during fetching, see if it's the
        position we want, when we find, we do a fetch() on it

        :param pk: The primary key from the position we want.
        :type pk: str

        :returns: A number pertaining to the kmer position.
        """

        return int(pk.split('|')[0])

    @staticmethod
    def get_lazy_ref_field_by_pos(pos_num: int, positions: list[LazyReferenceField]) -> LazyReferenceField:
        """
        A quick method to get the matching lazy ref field using position number.

        :param pos_num: The k-mer position number
        :param positions: A list of lazy ref fields of positions for a given job

        :type pos_num: int
        :type positions: List[LazyReferenceField]

        :returns: A LazyReferenceField for the provided k-mer position
        """

        for position_qs in positions:
            if not HelperMethods.get_id_from_pk(position_qs.pk) == pos_num:
                continue

            return position_qs

    @staticmethod
    def append_pos_to_pos_dict(pos_dict: dict, position: int) -> dict:
        """
        A quick method to get the matching lazy ref field using position number.

        :param pos_dict: The position dictionary after it was derived from the LazyReferenceField
        :param position: The position number.

        :type pos_dict: dict
        :type position: int

        :returns: A corrected dict containing the correct position number.
        """

        # Delete the pk which contains "posnum|job_id" because we do not need it
        del pos_dict['_id']

        # Add the position number
        pos_dict['position'] = position

        return pos_dict
