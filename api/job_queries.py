from flask_restful import abort
from mongoengine import DoesNotExist

from atat_single.models import Job, Result


class JobQueries(object):
    @staticmethod
    def get_job(jobid: str) -> Job:
        try:
            return Job.objects.get(_id=jobid)
        except DoesNotExist:
            abort(404, message=f'Job id {jobid} does not exist.')

    @staticmethod
    def get_result(jobid: str) -> Result:
        job = JobQueries.get_job(jobid)

        if job.status != 'FINISHED':
            abort(102, message=f'Job id {jobid} is {job.status}.')

        return job.results
