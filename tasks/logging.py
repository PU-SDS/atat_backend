from datetime import datetime


class Logging(object):
    @staticmethod
    def make_log_entry(context: str, msg: str) -> dict:
        """
            :param context: The heading of the log entry.
            :param msg: The actual message to include in the entry.

            :type context: str
            :type msg: str

            :returns: A dictionary  pertaining to the log entry.
        """

        timestamp = datetime.now().strftime("%Y/%M/%d at %H-%M-%S-%f")

        return {'context': context, 'msg': msg, 'timestamp': timestamp}
