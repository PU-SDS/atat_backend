from celery import Celery
from .settings import RabbitMQ, MongoDB

task_routes = {
    'ATAT': {'queue': 'ATAT'},
    'Hunana': {'queue': 'Hunana'},
    'Warehousing': {'queue': 'Warehousing'},
    'Job': {'queue': 'Job'}
}


app = Celery(
    broker=f'amqp://{RabbitMQ.USERNAME}:{RabbitMQ.PASSWORD}@{RabbitMQ.HOST}:{RabbitMQ.PORT}',
    backend=f'mongodb://{RabbitMQ.BACKEND_USERNAME}:{RabbitMQ.BACKEND_PASSWORD}@{MongoDB.HOST}:{MongoDB.PORT}/'
            f'{RabbitMQ.BACKEND_SOURCE}?authSource={RabbitMQ.BACKEND_SOURCE}'
)

app.autodiscover_tasks(['atat_single.tasks'])
app.conf.task_routes = task_routes

if __name__ == "__main__":
    app.start()
