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
    JOBS_FOLDER = '/home/shant/jobs'
    SPLIT_SEQS_SUBFOLDER = 'split_sequences'


class FileNames(Config):
    ENTROPY_CHART = 'entropy.jpg'
    VIOLIN_CHART = 'violin.jpg'
    JSON_POSITIONS = 'positions.json'
    JSON_ATAT_VARIANTS = 'variants.json'


class Containers(Config):
    MONGODB = 'mongodb'


class MongoDBs(Config):
    MONGO_POSITIONS_DB = 'positions'
    MONGO_VARIANTS_DB = 'variants'
    MONGO_AUTH_DB = 'f'


class FlaskConfig(Config):
    UPLOADS_DEFAULT_DEST = LocalPaths.JOBS_FOLDER
