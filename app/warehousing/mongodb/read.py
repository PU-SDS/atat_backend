import pandas as pd

from pymongo  import MongoClient
from app.warehousing.mongodb import MongoConstants

# Imported for type annotation
from pandas import DataFrame


class MongoDBRead(object):
    def __init__(self, job_id):
        mongo_client = MongoClient(host=MongoConstants.MONGO_HOST,
                                        port=MongoConstants.MONGO_PORT,
                                        username=MongoConstants.MONGODB_USR,
                                        password=MongoConstants.MONGODB_PWD,
                                        authSource=MongoConstants.MONGO_AUTH_DB,
                                        )[MongoConstants.MONGO_POSITIONS_DB]

        self.host_collection = mongo_client[f'{job_id}.source']
        self.reservoir_collection = mongo_client[f'{job_id}.reservoir']

    def get_all_positions(self, columns: dict) -> DataFrame:
        pd.options.plotting.backend = "plotly"

        positions_host = self.host_collection.find({}, columns)
        positions_reservoir = self.reservoir_collection.find({}, columns)

        df_host = pd.DataFrame(list(positions_host))
        df_reservoir = pd.DataFrame(list(positions_reservoir))

        merged_df = pd.merge(df_host, df_reservoir, on='position', suffixes=('-host', '-reservoir'))

        return merged_df

    def get_variants(self, position):
        aggregate_pipeline = [{'$match': {'position': position}},
            {'$unwind': '$sequences'},
            {'$project': {
                'sequences.position': 1,
                'sequences.sequence': 1,
                'sequences.count': 1,
                'sequences.conservation': 1,
                'sequences.motif_long': 1,
                'sequences.id': 1,
                'sequences.strain': 1,
                'sequences.country': 1,
                'sequences.host': 1,
                '_id': 0
            }}]

        variants_host = self.host_collection.aggregate(aggregate_pipeline)
        variants_reservoir = self.host_collection.aggregate(aggregate_pipeline)

        return [x['sequences'] for x in variants_host], [x['sequences'] for x in variants_reservoir]
