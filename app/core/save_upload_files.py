from app.core import GenerateJobId, Constants
from flask import current_app as app
from os.path import join
from os import makedirs


class SaveUploads(object):
    def __init__(self, host_file, reservoir_file):
        """
            Creates the required folders for the job and generates a unique job id

            Args:
                host_file: A list of dictionaries for each position of kmers in the host sequences
                reservoir_file: A list of dictionaries for each position of kmers in the reservoir sequences
        """

        self.host_file = host_file
        self.reservoir_file = reservoir_file

    def save(self) -> str:
        """
            Creates the required folder for the job and returns the job id

            Returns: Job id
        """
        if not self._is_filetype_valid(self.host_file.filename) or not \
                self._is_filetype_valid(self.host_file.filename):
            return False

        # Generate a unique job id
        job_id = GenerateJobId().generate()

        # Create our job path, folder and results sub folder
        job_folder = join(app.config['JOBS_FOLDER'], job_id)
        makedirs(join(job_folder, app.config['RESULTS_SUB_FOLDER']))

        # Save the files inside a unique folder inside the default job folder
        self.host_file.save(join(job_folder, Constants.HOST_SEQUENCE_FILENAME))
        self.reservoir_file.save(join(job_folder, Constants.RESERVOIR_SEQUENCE_FILENAME))

        # Return the job id
        return job_id

    def save2(self) -> str:
        """
            Creates the required folder for the job and returns the job id

            Returns: Job id
        """

        # Generate a unique job id
        job_id = GenerateJobId().generate()


        # Create our job path, folder and results sub folder
        job_folder = join(app.config['JOBS_FOLDER'], job_id)
        makedirs(join(job_folder, app.config['RESULTS_SUB_FOLDER']))

        with open(join(job_folder, Constants.HOST_SEQUENCE_FILENAME), 'wb') as host_file, \
                open(join(job_folder, Constants.RESERVOIR_SEQUENCE_FILENAME), 'wb') as reservoir_file:
            host_file.write(self.host_file)
            reservoir_file.write(self.reservoir_file)

        # Return the job id
        return job_id

    @classmethod
    def _is_filetype_valid(cls, filename):
        try:
            extension = filename.rsplit('.', 1)[1].lower()
        except:
            # TODO: Throw custom error
            pass

        if extension in Constants.ALLOWED_FILES:
            return True

        return False
