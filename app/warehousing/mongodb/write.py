from pymongo import MongoClient
from app.warehousing.mongodb.constants import MongoConstants


class MongoDBWrite(object):

    def __init__(self, job_id, source, reservoir):
        """
            Defines methods for saving JSON data and Chart images to disk

            :param job_id: The ID of the current job
            :type job_id: str
        """

        self.job_id = job_id
        self.source = source
        self.reservoir = reservoir
        self.mongo_client = MongoClient(host='localhost',
                                        port=27017,
                                        username=MongoConstants.MONGODB_USR,
                                        password=MongoConstants.MONGODB_PWD,
                                        authSource=MongoConstants.MONGO_AUTH_DB
                                        )

    def save_positions(self):
        """
            Saves the JSON encoded position objects to the database

            :return: A tuple containing the return code and output
        """

        positions_db = self.mongo_client[MongoConstants.MONGO_POSITIONS_DB]
        source_collection = positions_db[f'{self.job_id}.source']
        reservoir_collection = positions_db[f'{self.job_id}.reservoir']

        source_status = source_collection.insert_many(self.source)
        reservoir_status = reservoir_collection.insert_many(self.reservoir)

        return source_status, reservoir_status
