class JobExists(Exception):
    def __init__(self):
        """
        Raised when the job id provided by the user at the viva create job endpoint already exists in the database.
        """

        super(JobExists, self).__init__("Job already exists")
