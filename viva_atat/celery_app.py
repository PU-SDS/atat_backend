from celery import Celery
from .settings import ResourceSettings

settings = ResourceSettings()
task_routes = {
    'ATAT': {'queue': 'ATAT'},
    'Hunana': {'queue': 'Hunana'},
    'Warehousing': {'queue': 'Warehousing'},
    'Job': {'queue': 'Job'},
}

app = Celery(
    broker=f'amqp://{settings.rabbitmq_username}:{settings.rabbitmq_password}@{settings.rabbitmq_host}:'
    f'{settings.rabbitmq_port}',
    backend=f'mongodb://{settings.mongo_rabbitmq_username}:{settings.mongo_rabbitmq_password}@{settings.mongo_host}:'
    f'{settings.mongo_port}/{settings.mongo_rabbitmq_database}?replicaSet=rs0&authSource='
    f'{settings.mongo_rabbitmq_database}',
)

app.autodiscover_tasks(['viva_atat.tasks'])
app.conf.task_routes = task_routes

if __name__ == "__main__":
    app.start()
