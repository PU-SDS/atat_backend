import hashlib
from datetime import datetime

from mongoengine import DoesNotExist

from ..models import Job


class Logging(object):
    @staticmethod
    def make_log_entry(jobid: str, context: str, msg: str, status: str = None):
        """
            :param context: The heading of the log entry.
            :param msg: The actual message to include in the entry.
            :param status: The status of the job.

            :type context: str
            :type msg: str
            :type status: str
        """

        timestamp = datetime.now().strftime("%Y/%M/%d at %H-%M-%S-%f")
        hashx = hashlib.md5('-'.join((msg, timestamp)).encode('utf-8')).hexdigest()
        entry = {'hash': hashx, 'context': context, 'msg': msg, 'timestamp': timestamp}

        try:
            job = Job.objects.get(_id=jobid)
        except DoesNotExist:
            Job(_id=jobid, log=[entry], status=status).save()
            return

        if status:
            job.status = status

        job.log.append(entry)
        job.save()
