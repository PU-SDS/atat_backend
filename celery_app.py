from celery import Celery

task_routes = {
    'ATAT': {'queue': 'ATAT'},
    'Hunana': {'queue': 'Hunana'},
    'Warehousing': {'queue': 'Warehousing'},
    'Job': {'queue': 'Job'}
}

app = Celery(
    broker='amqp://defaultrabbit:defaultrabbitpass@localhost:5672',
    backend='mongodb://rabbitmongousr:rabbitmongopwd@127.0.0.1:27017/rabbit_backend?authSource=rabbit_backend'
)

app.autodiscover_tasks(['atat_single.tasks'])
app.conf.task_routes = task_routes
# app.conf.task_default_queue = 'Warehousing'

if __name__ == "__main__":
    app.start()
