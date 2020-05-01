from celery import Celery

app = Celery('flua2h',backend='rpc://',broker='pyamqp://rabbitmq:rabbitmq@viva-rabbitmq//',include=['flua2h.tasks'])
#app = Celery('flua2h',backend='rpc://',broker='pyamqp://rabbitmq:rabbitmq@viva-rabbitmq//',include=['flua2h.tasks','common.jobmanagement','common.sequence','common.hunana'])
