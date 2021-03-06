from os import environ
from mongoengine import connect, StringField, IntField, ListField, EmbeddedDocumentListField, PULL, \
    EmbeddedDocumentField, SequenceField, FloatField
from mongoengine.base import EmbeddedDocumentList
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


class Job(Document):
    """
        This is the model for a basic ATAT job.
    """

    _id = StringField(primary_key=True, required=True)
    log = ListField(required=True)
    results = FollowReferenceField(Result, required=False)
