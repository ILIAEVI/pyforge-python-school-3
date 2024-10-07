from celery import Celery


celery = Celery(
    'tasks',
    broker='redis://redis:6379/0',
    backend='redis://redis:6379/0',
)

celery.conf.update(
    broker_connection_retry_on_startup=True
)

celery.autodiscover_tasks(['src.tasks'])