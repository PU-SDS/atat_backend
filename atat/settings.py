from os import environ


class RabbitMQ(object):
    class Dev(object):
        HOST = environ.get('RABBITMQ_HOST')
        PORT = environ.get('RABBITMQ_PORT')

        USERNAME = environ.get('RABBITMQ_USERNAME')
        PASSWORD = environ.get('RABBITMQ_PASSWORD')

        BACKEND_SOURCE = environ.get('RABBITMQ_BACKEND_SOURCE')
        BACKEND_USERNAME = environ.get('RABBITMQ_BACKEND_USERNAME')
        BACKEND_PASSWORD = environ.get('RABBITMQ_BACKEND_PASSWORD')

    class Prod(object):
        HOST = environ.get('RABBITMQ_ATAT_SERVICE_HOST')
        PORT = 5672

        USERNAME = environ.get('username')
        PASSWORD = environ.get('password')

        BACKEND_ = environ.get('RABBITMQ_BACKEND_SOURCE')
        BACKEND_USERNAME = environ.get('RABBITMQ_BACKEND_USERNAME')
        BACKEND_PASSWORD = environ.get('RABBITMQ_BACKEND_PASSWORD')


class MongoDB(object):
    HOST = environ.get('MONGODB_HOST')
    PORT = environ.get('MONGODB_PORT')

    USERNAME = environ.get('MONGODB_USERNAME')
    PASSWORD = environ.get('MONGODB_PASSWORD')

    PROJECT_DATABASE = environ.get('MONGO_PROJECT_DATABASE')
