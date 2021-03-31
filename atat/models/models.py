from mongoengine import connect, StringField, IntField, ListField, EmbeddedDocumentListField, FloatField, QuerySet
from mongoengine_goodjson import Document, EmbeddedDocument, FollowReferenceField

from ..settings import MongoDB

connect(MongoDB.PROJECT_DATABASE,
        authentication_source=MongoDB.PROJECT_DATABASE,
        username=MongoDB.USERNAME,
        password=MongoDB.PASSWORD,
        host=MongoDB.HOST,
        port=MongoDB.PORT
        )


class HunanaPosition(EmbeddedDocument):
    """
        This is the model for a Hunana result position. For both source, and reservoir.
    """

    position = IntField(required=True)
    entropy = FloatField(required=True)
    supports = IntField(required=True)
    variants = ListField(required=True)


class Switch(EmbeddedDocument):
    """
        This is the model for a motif switch.
    """

    class Motifs(object):
        """
            Defines various motif names.
        """

        MOTIFS = {
            'INDEX': 'Index',
            'MAJOR': 'Major',
            'MINOR': 'Minor',
            'UNIQUE': 'Unique'
        }

    position = IntField(required=True)
    sequence = StringField(required=True)
    fromx = StringField(required=True, choices=Motifs.MOTIFS.keys())
    to = StringField(required=True, choices=Motifs.MOTIFS.keys())


class Result(Document):
    """
        This is the model for a ATAT job result.
    """

    source = EmbeddedDocumentListField(HunanaPosition)
    reservoir = EmbeddedDocumentListField(HunanaPosition)
    switches = EmbeddedDocumentListField(Switch)

    class ResultQuerySet(QuerySet):
        def get_grouped_position(self, idx: str, position: int):
            result = self.get(id=idx)

            source_position = result.source.get(position=position).to_mongo().to_dict()
            reservoir_position = result.reservoir.get(position=position).to_mongo().to_dict()

            return {"source": source_position, "reservoir": reservoir_position}

    meta = {'queryset_class': ResultQuerySet}


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
