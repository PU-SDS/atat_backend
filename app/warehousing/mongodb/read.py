from typing import Tuple, List

import pandas as pd

from pymongo import MongoClient
from bson.json_util import dumps
from app.warehousing.mongodb.constants import MongoConstants, JobStatuses

# Imported for type annotation
from pandas import DataFrame


class MongoDBRead(object):
    def __init__(self, job_id):
        self.mongo_client = MongoClient(host=MongoConstants.MONGO_HOST,
                                   port=MongoConstants.MONGO_PORT,
                                   username=MongoConstants.MONGODB_USR,
                                   password=MongoConstants.MONGODB_PWD,
                                   authSource=MongoConstants.MONGO_AUTH_DB,
                                   )

        self.job_id = job_id

    def get_all_positions(self, columns: dict) -> DataFrame:
        pd.options.plotting.backend = "plotly"

        positions_db = self.mongo_client[MongoConstants.MONGO_POSITIONS_DB]

        host_collection = positions_db[f'{self.job_id}.source']
        reservoir_collection = positions_db[f'{self.job_id}.reservoir']

        positions_host = host_collection.find({}, columns)
        positions_reservoir = reservoir_collection.find({}, columns)

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

        positions_db = self.mongo_client[MongoConstants.MONGO_POSITIONS_DB]

        host_collection = positions_db[f'{self.job_id}.source']
        reservoir_collection = positions_db[f'{self.job_id}.reservoir']

        variants_host = host_collection.aggregate(aggregate_pipeline)
        variants_reservoir = reservoir_collection.aggregate(aggregate_pipeline)

        return [x['sequences'] for x in variants_host], [x['sequences'] for x in variants_reservoir]

    def get_max_positions(self):
        positions_db = self.mongo_client[MongoConstants.MONGO_POSITIONS_DB]

        host_collection = positions_db[f'{self.job_id}.source']

        max_position = list(host_collection.aggregate(
            [
                {'$group': {
                    '_id': 0,
                    'max': {'$max': '$position'}
                }}
            ]
        ))[0].get('max')

        return max_position

    def get_position_statistics(self, position: int) -> Tuple:
        positions_db = self.mongo_client[MongoConstants.MONGO_POSITIONS_DB]

        host_collection = positions_db[f'{self.job_id}.source']
        reservoir_collection = positions_db[f'{self.job_id}.reservoir']

        host_statistics = host_collection.find({'position': position}, {'entropy': 1, 'supports': 1,
                                                                             'variants': 1, '_id': 0})

        reservoir_statistics = reservoir_collection.find({'position': position}, {'entropy': 1, 'supports': 1,
                                                                                       'variants': 1, '_id': 0})

        return list(host_statistics), list(reservoir_statistics)

    def get_job_status(self):
        positions_db = self.mongo_client[MongoConstants.MONGO_JOB_STATUS_DB]

        job_status_collection = positions_db[self.job_id]
        status = list(job_status_collection.find({}, {'status': 1}))

        if not status:
            return JobStatuses.JOB_NONEXIST

        return status[0].get('status')

    def get_motifs(self, source_motif: str, reservoir_motif: str) -> Tuple[List, List]:
        positions_db = self.mongo_client[MongoConstants.MONGO_POSITIONS_DB]

        host_collection = positions_db[f'{self.job_id}.source']
        reservoir_collection = positions_db[f'{self.job_id}.reservoir']

        # Get motifs from source matching the provided source motif
        source_motifs = host_collection.aggregate(
            self._get_motifs_pipeline(source_motif)
        )

        # Get motifs from reservoir matching the provided reservoir motif
        reservoir_motifs = reservoir_collection.aggregate(
            self._get_motifs_pipeline(reservoir_motif)
        )

        return list(source_motifs), list(reservoir_motifs)

    def get_json(self) -> Tuple:
        positions_db = self.mongo_client[MongoConstants.MONGO_POSITIONS_DB]

        host_collection = positions_db[f'{self.job_id}.source']
        reservoir_collection = positions_db[f'{self.job_id}.reservoir']

        host_results = host_collection.find({}, {'_id': 0, 'variants_flattened': 0})
        reservoir_results = reservoir_collection.find({}, {'_id': 0, 'variants_flattened': 0})

        return host_results, reservoir_results

    @classmethod
    def _get_motifs_pipeline(cls, motif_short: str) -> list:
        # A centralized place to derive the pipeline for motif searching
        pipeline = [
                {'$unwind': '$sequences'},
                {'$match': {'sequences.motif_short': motif_short}},
                {'$group': {'_id': '$sequences.position', 'variants': {'$push': '$$ROOT'}}},
                {'$project': {'variants.sequences': 1}}
            ]

        return pipeline
