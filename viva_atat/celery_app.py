from celery import Celery
from .settings import ResourceSettings

settings = ResourceSettings()

task_routes = {
    'atat': {'queue': 'atat'},
    'dima': {'queue': 'dima'},
    'warehousing': {'queue': 'warehousing'},
    'job': {'queue': 'job'},
}

app = Celery(
    broker=f'amqp://{settings.rabbitmq_username}:'
    f'{settings.rabbitmq_password}@'
    f'{settings.rabbitmq_host}:'
    f'{settings.rabbitmq_port}',
    backend=f'mongodb://{settings.mongo_rabbitmq_username}:'
    f'{settings.mongo_rabbitmq_password}@'
    f'{settings.mongo_host}:'
    f'{settings.mongo_port}/'
    f'{settings.mongo_rabbitmq_database}'
    f'?authSource={settings.mongo_rabbitmq_database}'
    + (f'&replicaSet={settings.mongo_replicaset_name}' if settings.env_state == "prod" else ''),
)

app.autodiscover_tasks(['viva_atat.core.tasks'])
app.conf.task_routes = task_routes

if __name__ == "__main__":
    app.start()
