from celery import Celery

app = Celery(
    broker='amqp://defaultrabbit:defaultrabbitpass@localhost:5672',
    backend='mongodb://rabbitmongousr:rabbitmongopwd@localhost:27017/rabbit_backend?authSource=rabbit_backend'
)

app.autodiscover_tasks(['atat_single.tasks'])

if __name__ == "__main__":
    app.start()
