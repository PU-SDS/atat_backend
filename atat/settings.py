from os import environ


class MongoDB(object):
    class Dev(object):
        HOST = 'mongodb-svc.atat'
        PORT = 27017

        USERNAME = 'atat'
        PASSWORD = environ.get('MONGODB_ATAT_PASSWORD')

        PROJECT_DATABASE = 'atat'

    class Prod(object):
        HOST = 'mongodb-headless.default.svc.cluster.local'
        PORT = 27017

        USERNAME = 'atat'
        PASSWORD = environ.get('MONGODB_ATAT_PASSWORD')

        PROJECT_DATABASE = 'atat'

    @staticmethod
    def get_connection_string(env: str, mongo_replicas: int = 3, podname_prefix: str = 'mongodb-'):
        host = MongoDB.Prod.HOST if env == 'prod' else MongoDB.Dev.HOST
        port = MongoDB.Prod.PORT if env == 'prod' else MongoDB.Dev.PORT

        username = MongoDB.Prod.USERNAME if env == 'prod' else MongoDB.Dev.USERNAME
        password = MongoDB.Prod.PASSWORD if env == 'prod' else MongoDB.Dev.PASSWORD
        database = MongoDB.Prod.PROJECT_DATABASE if env == 'prod' else MongoDB.Dev.PROJECT_DATABASE

        con_strings = []

        for replica in range(mongo_replicas - 1):
            con_string = f'{username}:{password}@{podname_prefix}{replica}.{host}:{port}/{database}/?authSource={database}'
            con_strings.append(con_string)

        return f'mongodb://{",".join(con_strings)}/?replicaSet=mongodb&ssl=false&tls=false'


class RabbitMQ(object):
    class Dev(object):
        HOST = environ.get('RABBITMQ_ATAT_SERVICE_HOST')
        PORT = environ.get('RABBITMQ_ATAT_SERVICE_PORT')

        USERNAME = environ.get('username')
        PASSWORD = environ.get('password')

        BACKEND_DATABASE = 'rabbitmq'
        BACKEND_USERNAME = 'rabbitmq'
        BACKEND_PASSWORD = environ.get('MONGODB_RABBITMQ_PASSWORD')

    class Prod(object):
        HOST = 'rabbitmq.atat-dev.svc'
        PORT = 5672

        USERNAME = 'rabbitmq'
        PASSWORD = environ.get('RABBITMQ_BROKER_PASSWORD')

        BACKEND_DATABASE = 'rabbitmq'
        BACKEND_USERNAME = 'rabbitmq'
        BACKEND_PASSWORD = environ.get('MONGODB_RABBITMQ_PASSWORD')

    @staticmethod
    def get_backend_connection_string(env: str, mongo_replicas: int = 3, podname_prefix: str = 'mongodb-'):
        host = MongoDB.Prod.HOST if env == 'prod' else MongoDB.Dev.HOST
        port = MongoDB.Prod.PORT if env == 'prod' else MongoDB.Dev.PORT

        username = RabbitMQ.Prod.BACKEND_USERNAME if env == 'prod' else RabbitMQ.Dev.BACKEND_USERNAME
        password = RabbitMQ.Prod.BACKEND_PASSWORD if env == 'prod' else RabbitMQ.Dev.BACKEND_PASSWORD
        database = RabbitMQ.Prod.BACKEND_DATABASE if env == 'prod' else RabbitMQ.Dev.BACKEND_DATABASE

        con_strings = []

        for replica in range(mongo_replicas - 1):
            con_string = f'{username}:{password}@{podname_prefix}{replica}.{host}:{port}/{database}'
            con_strings.append(con_string)

        return f'mongodb://{",".join(con_strings)}/?replicaSet=mongodb&ssl=false&tls=false&authSource={database}'
