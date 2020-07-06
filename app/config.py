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
    MAAFT = 'maaft'
    TRIMAL = 'trimal'


class LocalPaths(Config):
    JOBS_FOLDER = '/home/shant/jobs'
    SPLIT_SEQS_SUBFOLDER = 'split_sequences'
