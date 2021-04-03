from celery import Celery
from .settings import RabbitMQ, MongoDB

task_routes = {
    'ATAT': {'queue': 'ATAT'},
    'Hunana': {'queue': 'Hunana'},
    'Warehousing': {'queue': 'Warehousing'},
    'Job': {'queue': 'Job'}
}


app = Celery(
    broker=f'amqp://guest:guest@{RabbitMQ.Prod.HOST}:{RabbitMQ.Prod.PORT}'
)

app.autodiscover_tasks(['atat.tasks'])
app.conf.task_routes = task_routes

if __name__ == "__main__":
    app.start()
