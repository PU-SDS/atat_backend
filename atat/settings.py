from os import environ


class RabbitMQ(object):
    HOST = environ.get('RABBITMQ_HOST')
    PORT = environ.get('RABBITMQ_PORT')

    USERNAME = environ.get('RABBITMQ_USERNAME')
    PASSWORD = environ.get('RABBITMQ_PASSWORD')

    BACKEND_SOURCE = environ.get('RABBITMQ_BACKEND_SOURCE')
    BACKEND_USERNAME = environ.get('RABBITMQ_BACKEND_USERNAME')
    BACKEND_PASSWORD = environ.get('RABBITMQ_BACKEND_PASSWORD')


class MongoDB(object):
    HOST = environ.get('MONGODB_HOST')
    PORT = environ.get('MONGODB_PORT')

    USERNAME = environ.get('MONGODB_USERNAME')
    PASSWORD = environ.get('MONGODB_PASSWORD')

    PROJECT_DATABASE = environ.get('MONGO_PROJECT_DATABASE')
