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
    JOBS_FOLDER = 'C:\Jobs\\'
    SPLIT_SEQS_SUBFOLDER = 'split_sequences'
