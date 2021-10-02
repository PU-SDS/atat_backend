from celery import Celery
from .settings import ResourceSettings

settings = ResourceSettings()
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

if __name__ == "__main__":
    app.start()
