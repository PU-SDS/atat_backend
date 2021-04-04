from celery import Celery
from .settings import RabbitMQ, MongoDB

task_routes = {
    'ATAT': {'queue': 'ATAT'},
    'Hunana': {'queue': 'Hunana'},
    'Warehousing': {'queue': 'Warehousing'},
    'Job': {'queue': 'Job'}
}

app = Celery(
    broker=f'amqp://{RabbitMQ.Prod.USERNAME}:{RabbitMQ.Prod.PASSWORD}@{RabbitMQ.Prod.HOST}:{RabbitMQ.Prod.PORT}',
    backend=f'mongodb://{RabbitMQ.Prod.BACKEND_USERNAME}:{RabbitMQ.Prod.BACKEND_PASSWORD}@{MongoDB.Prod.HOST}:'
            f'{MongoDB.Prod.PORT}/{RabbitMQ.Prod.BACKEND_DATABASE}?replicaSet=mongodb&authSource='
            f'{RabbitMQ.Prod.BACKEND_DATABASE}'
)

app.autodiscover_tasks(['atat.tasks'])
app.conf.task_routes = task_routes

if __name__ == "__main__":
    app.start()
