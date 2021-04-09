import hashlib
from datetime import datetime

from .constants import JOB_ID_GLOBAL
from ..models import Job


class Logging(object):
    @staticmethod
    def make_log_entry(context: str, msg: str):
        """
            :param context: The heading of the log entry.
            :param msg: The actual message to include in the entry.

            :type context: str
            :type msg: str
        """

        job = Job.objects.get(_id=JOB_ID_GLOBAL)

        timestamp = datetime.now().strftime("%Y/%M/%d at %H-%M-%S-%f")
        hashx = hashlib.md5('-'.join((msg, timestamp)).encode('utf-8')).hexdigest()
        entry = {'hash': hashx, 'context': context, 'msg': msg, 'timestamp': timestamp}

        job.log.append(entry)
