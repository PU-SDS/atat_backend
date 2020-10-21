from os import path, mkdir
from flask import current_app as app
import jsonpickle


class FlatFile(object):
    def __init__(self, job_id: str):
        """
            Defines methods for saving JSON data and Chart images to disk

            :param job_id: The ID of the current job
            :type job_id: str
        """

        self.job_path = path.join(app.config['JOBS_FOLDER'], job_id)

    def create_directory(self) -> str:
        """
            Create the job folder

            :return: The full path to the job folder
        """

        try:
            mkdir(self.job_path)
        except OSError:
            # TODO: Custom Error Handler here
            print(f'Could not create job folder for ID {self.job_path}')

        return self.job_path

    def save_positions_json(self, data: list) -> str:
        """
            Saves a given list of Position objects to the disk as JSON

            :param data: A list of Position objects
            :type data: list

            :return: The full path to the saved JSON file
        """

        save_path = path.join(self.job_path, app.config['JSON_POSITIONS'])

        self._save_json(save_path, data)

        return save_path

    def save_atat_json(self, data: list) -> str:
        """
            Saves a given list of ATATData objects to the disk as JSON

            :param data: A list of ATATData objects
            :type data: list

            :return: The full path to the saved JSON file
        """

        save_path = path.join(self.job_path, app.config['JSON_ATAT_VARIANTS'])

        self._save_json(save_path, data)

        return save_path

    def _save_json(self, save_path: str, data: list):
        """
            Saves the JSON file to the disk

            :param save_path: The full path to save the file at
            :param data: A list containing Position or ATATData objects

            :type save_path: str
            :type data: list
        """

        json_data = jsonpickle.encode(data, unpicklable=False)

        with open(save_path, 'w') as f:
            f.write(json_data)
