from flask import current_app as app
from os import environ


class MongoConstants(object):
    MONGODB_USR = environ['MONGO_USER']
    MONGODB_PWD = environ['MONGO_PWD']

    MONGO_AUTH_DB = app.config['MONGO_AUTH_DB']
    MONGO_HOST = app.config['MONGO_HOST']
    MONGO_PORT = app.config['MONGO_PORT']

    MONGO_POSITIONS_DB = app.config['MONGO_POSITIONS_DB']
    MONGO_JOB_STATUS_DB = app.config['MONGO_JOB_STATUS_DB']


class JobStatuses(object):
    JOB_PROCESSING = 'Processing'
    JOB_COMPLETED = 'Completed'
    JOB_FAILED = 'Failed'
    JOB_NONEXIST = 'Non-Existant'
