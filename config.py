from flask import Config


class ProductionConfig(Config):
    DEBUG = False
    TESTING = False


class DevelopmentConfig(Config):
    DEBUG = True
    TESTING = True


class TestingConfig(Config):
    TESTING = True


class BioApps(Config):
    MAFFT = 'mafft'
    TRIMAL = 'trimal'


class LocalPaths(Config):
    JOBS_FOLDER = 'C:\\Jobs\\'
    SPLIT_SEQS_SUBFOLDER = 'split_sequences'
    RESULTS_SUB_FOLDER = 'results'


class FileNames(Config):
    ENTROPY_CHART = 'entropy.jpg'
    VIOLIN_CHART = 'violin.jpg'
    JSON_POSITIONS = 'positions.json'
    JSON_ATAT_VARIANTS = 'variants.json'


class Containers(Config):
    MONGODB = 'mongodb'


class MongoDBs(Config):
    MONGO_POSITIONS_DB = 'positions'
    MONGO_JOB_STATUS_DB = 'job_status'
    MONGO_VARIANTS_DB = 'variants'
    MONGO_AUTH_DB = 'positions'
    MONGO_HOST = 'localhost'
    MONGO_PORT = 27017


class CeleryConfig(Config):
    CELERY_BROKER_URL = 'pyamqp://guest@localhost//'