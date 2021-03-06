from os import environ
from mongoengine import connect, StringField, IntField, ListField, EmbeddedDocumentListField, FloatField, QuerySet
from mongoengine_goodjson import Document, EmbeddedDocument, FollowReferenceField

DATABASE_NAME = environ.get('MONGO_INITDB_DATABASE')
USERNAME = environ.get('MONGO_ATAT_USERNAME')
PASSWORD = environ.get('MONGO_ATAT_PASSWORD')

connect(DATABASE_NAME,
        authentication_source=DATABASE_NAME,
        username=USERNAME,
        password=PASSWORD
        )


class HunanaPosition(EmbeddedDocument):
    """
        This is the model for a Hunana result position. For both source, and reservoir.
    """

    position = IntField(required=True)
    entropy = FloatField(required=True)
    supports = IntField(required=True)
    variants = ListField(required=True)
    kmer_types = ListField(required=True)


class Result(Document):
    """
        This is the model for a ATAT job result.
    """

    source = EmbeddedDocumentListField(HunanaPosition)
    reservoir = EmbeddedDocumentListField(HunanaPosition)

    class ResultQuerySet(QuerySet):
        def get_grouped_position(self, idx: str, position: int):
            result = self.get(id=idx)

            source_position = result.source.get(position=position).to_mongo().to_dict()
            reservoir_position = result.reservoir.get(position=position).to_mongo().to_dict()

            return {"source": source_position, "reservoir": reservoir_position}

    meta = {'queryset_class':  ResultQuerySet}


class Job(Document):
    """
        This is the model for a basic ATAT job.
    """

    statuses = {
        'STARTED': 'STARTED',
        'RUNNING': 'RUNNING',
        'FAILED': 'FAILED',
        'FINISHED': 'FINISHED'
    }

    @property
    def get_status(self):
        return self.statuses.get(self.status)

    _id = StringField(primary_key=True, required=True)
    status = StringField(required=True, choices=statuses.keys(), max_length=8)
    log = ListField(required=True)
    results = FollowReferenceField(Result, required=False)
