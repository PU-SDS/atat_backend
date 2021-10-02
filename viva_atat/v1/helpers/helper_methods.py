from ...core.models import JobDBModel, Parameters
from ..models import CreateJobParameters


class HelperMethods(object):
    @staticmethod
    def register_job(parameters: CreateJobParameters) -> str:
        """
        Given a valid set of parameters originating from a valid job create request, registers it as a new job on the
        database and returns the auto-generated job id.

        :param parameters: The parameters originating from a valid request to the job creation endpoint.
        :type parameters: CreateJobParameters

        :returns: The auto-generated job id.
        """

        parameters = Parameters(kmer_length=parameters.kmer_length, header_format=parameters.header_format).save()
        job = JobDBModel(parameters=parameters).save()

        return job.id

    @staticmethod
    def get_job(job_id: str) -> JobDBModel:
        """
        Given a valid job id, returns an instance of the job model pertaining to the job id.

        :param job_id: The ID of the job to retrieve.
        :type job_id: str

        :returns: An instance of a job model pertaining to the job id.
        """

        return JobDBModel.objects.get(id=job_id)
