from datetime import datetime
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
    UUIDField,
    DateTimeField,
    DictField,
    Document,
    EmbeddedDocument,
    ReferenceField,
)

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
    """
    Possible statuses for ATAT jobs.
    """

    CREATED: str = 'Created'
    RUNNING: str = 'Running'
    FAILED: str = 'Failed'
    COMPLETED: str = 'Completed'


class MotifClasses(Enum):
    """
    Motif classes as classified by DiMA.
    """

    index: str = 'Index'
    major: str = 'Major'
    minor: str = 'Minor'
    unique: str = 'Unique'


class LogMessageFlags(Enum):
    """
    The possible flags for a log message to indicate severity.
    """

    INFO: str = 'Information'
    ERROR: str = 'Error'
    WARNING: str = 'Warning'
    SUCCESS: str = 'Success'


class LogMessages(Enum):
    JOB_CREATED = 'Job created.'
    JOB_RUNNING = 'Job is running.'
    RUN_MOTIF_CLASSIFICATION = 'Running motif classification using DiMA.'
    MOTIF_CLASSIFICATION_ERROR = 'Error occurred while running DiMA.'
    MOTIF_CLASSIFICATION_COMPLETE = 'Completed motif classification.'
    RUN_TRANSMISSIBILITY_ANALYSIS = 'Running transmissibility analysis.'
    TRANSMISSIBILITY_ANALYSIS_ERROR = 'Error occurred during transmissibility analysis.'
    SAVING_TO_DB = 'Saving results to database.'
    ADDED_MOTIF_RESULTS = 'Saved motif classification results on the database.'
    ADDED_MOTIF_SWITCH_RESULTS = 'Saved motif switch results on the database.'
    JOB_COMPLETED = 'Analysis completed.'
    JOB_FAILED = 'An error occurred. Job cancelled.'


class LogEntryDBModel(EmbeddedDocument):
    """
    The model for a log entry.
    """

    id = UUIDField(required=True, default=uuid4, binary=False, db_field='id')
    flag = EnumField(LogMessageFlags, required=True)
    timestamp = DateTimeField(required=True, default=datetime.now)
    message = EnumField(LogMessages, required=True)


class DimaVariant(EmbeddedDocument):
    """
    This is the model for a kmer variant.
    """

    sequence = StringField(required=True)
    count = IntField(required=True)
    incidence = FloatField(required=True)
    motif = StringField(required=True, choices=[motif.value for motif in MotifClasses])
    metadata = DictField(required=True)


class DimaPosition(EmbeddedDocument):
    """
    This is the model for a Hunana result position. For both host, and reservoir.
    """

    position = IntField(required=True)
    supports = IntField(required=True)
    variants = EmbeddedDocumentListField(DimaVariant, default=[])


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

    host = EmbeddedDocumentListField(DimaPosition, required=False)
    reservoir = EmbeddedDocumentListField(DimaPosition, required=False)
    switches = EmbeddedDocumentListField(Switch, required=False)

    class ResultQuerySet(QuerySet):
        def get_grouped_position(self, position: int):
            result = self.get()

            host_position = result.host.get(position=position).to_mongo().to_dict()
            reservoir_position = result.reservoir.get(position=position).to_mongo().to_dict()

            return {"host": host_position, "reservoir": reservoir_position}

    meta = {'queryset_class': ResultQuerySet}


class Parameters(Document):
    """
    The model for the job and DiMA parameters
    """

    kmer_length = IntField(required=True)
    header_format = ListField(required=True)
    email = StringField(required=False, default=None)


class JobDBModel(Document):
    """
    The model for an ATAT job.
    """

    class LoggerQuerySet(QuerySet):
        def update_log(self, flag: LogMessageFlags, msg: LogMessages):
            entry = LogEntryDBModel(flag=flag, message=msg)
            job = self.get()

            job.log.append(entry)
            job.save()

        def update_status(self, status: JobStatus):
            job = self.get()

            job.status = status
            job.save()

    id = StringField(required=True, default=lambda: str(uuid4()), primary_key=True)
    status = EnumField(JobStatus, default=JobStatus.CREATED)
    parameters = ReferenceField(Parameters, required=True)
    log = EmbeddedDocumentListField(LogEntryDBModel, required=False)
    results = ReferenceField(Results, required=False, default=Results().save())

    meta = {'collection': 'job', 'queryset_class': LoggerQuerySet}
