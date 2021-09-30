from enum import Enum
from uuid import uuid4

from mongoengine import (
    connect,
    StringField,
    IntField,
    ListField,
    EmbeddedDocumentListField,
    FloatField,
    QuerySet,
    EnumField,
)
from mongoengine_goodjson import Document, EmbeddedDocument, FollowReferenceField

from ...settings import ResourceSettings

settings = ResourceSettings()

connect(
    db=settings.mongo_atat_database,
    host=settings.mongo_host,
    port=settings.mongo_port,
    authentication_source=settings.mongo_atat_database,
    username=settings.mongo_atat_username,
    password=settings.mongo_atat_password,
    **{'replicaset': settings.mongo_replicaset_name} if settings.env_state == "prod" else {},
)


class JobStatus(Enum):
    pending: str = 'pending'
    running: str = 'running'
    failed: str = 'failed'
    completed: str = 'completed'


class MotifClasses(Enum):
    index: str = 'Index'
    major: str = 'Major'
    minor: str = 'Minor'
    unique: str = 'Unique'


class DimaPosition(EmbeddedDocument):
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

    position = IntField(required=True)
    sequence = StringField(required=True)
    source = EnumField(MotifClasses, required=True)
    target = EnumField(MotifClasses, required=True)


class Results(Document):
    """
    This is the model for a ATAT job result.
    """

    source = EmbeddedDocumentListField(DimaPosition)
    reservoir = EmbeddedDocumentListField(DimaPosition)
    switches = EmbeddedDocumentListField(Switch)

    class ResultQuerySet(QuerySet):
        def get_grouped_position(self, idx: str, position: int):
            result = self.get(id=idx)

            source_position = result.source.get(position=position).to_mongo().to_dict()
            reservoir_position = result.reservoir.get(position=position).to_mongo().to_dict()

            return {"source": source_position, "reservoir": reservoir_position}

    meta = {'queryset_class': ResultQuerySet}


class JobDBModel(Document):
    id = StringField(required=True, default=lambda: str(uuid4()), primary_key=True)
    status = EnumField(JobStatus, default=JobStatus.pending)
    log = ListField(required=False, null=True)
    results = FollowReferenceField(Results, required=False)

    meta = {'collection': 'job'}
