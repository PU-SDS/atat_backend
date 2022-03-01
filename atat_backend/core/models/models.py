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
    CASCADE,
    LazyReferenceField
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

        - Created
        - Running
        - Failed
        - Completed
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
    JOB_COMPLETED = 'Analysis completed.'
    JOB_FAILED = 'An error occurred. Job cancelled.'
    RUN_MOTIF_CLASSIFICATION = 'Running motif classification using DiMA.'
    MOTIF_CLASSIFICATION_ERROR = 'Error occurred while running DiMA.'
    MOTIF_CLASSIFICATION_COMPLETE = 'Completed motif classification.'
    RUN_TRANSMISSIBILITY_ANALYSIS = 'Running transmissibility analysis.'
    TRANSMISSIBILITY_ANALYSIS_ERROR = 'Error occurred during transmissibility analysis.'
    TRANSMISSIBILITY_ANALYSIS_COMPLETE = 'Transmissibility analysis completed.'
    SAVING_TO_DB = 'Saving results to database.'
    ADDED_MOTIF_RESULTS = 'Saved motif classification results in the database.'
    ADDED_MOTIF_SWITCH_RESULTS = 'Saved motif transmission results in the database.'


class LogEntryDBModel(EmbeddedDocument):
    """
    The model for a log entry.
    """

    id = UUIDField(required=True, default=uuid4, binary=False, db_field='id')
    flag = EnumField(LogMessageFlags, required=True)
    timestamp = DateTimeField(required=True, default=datetime.now)
    message = EnumField(LogMessages, required=True)


class DimaVariantDBModel(Document):
    """
    This is the model for a kmer variant.
    """

    sequence = StringField(required=True)
    count = IntField(required=True)
    incidence = FloatField(required=True)
    motif_long = EnumField(MotifClasses, required=True)
    metadata = DictField(required=True)

    meta = {'strict': False}


class DimaPositionDBModel(Document):
    """
    This is the model for a Hunana result position. For both dataset one, and dataset two.
    """

    _fields = None
    position = StringField(required=True, primary_key=True)
    support = IntField(required=True)
    variants = ListField(ReferenceField(DimaVariantDBModel, reverse_delete_rule=CASCADE))

    meta = {'strict': False}


class TransmissionDBModel(Document):
    """
    This is the model for a motif switch.
    """

    position = IntField(required=True)
    sequence = StringField(required=True)
    source = EnumField(MotifClasses, required=True)
    target = EnumField(MotifClasses, required=True)
    source_incidence = FloatField(required=True)
    target_incidence = FloatField(required=True)


class ResultsDBModel(Document):
    """
    This is the model for a ATAT job result.
    """

    dataset_one = ListField(LazyReferenceField(DimaPositionDBModel, reverse_delete_rule=CASCADE))
    dataset_two = ListField(LazyReferenceField(DimaPositionDBModel, reverse_delete_rule=CASCADE))
    switches = ListField(ReferenceField(TransmissionDBModel, reverse_delete_rule=CASCADE), default=[])


class ParametersDBModel(Document):
    """
    The model for the job and DiMA parameters
    """

    kmer_length = IntField(required=True)
    header_format = ListField(required=True)
    protein_name = StringField(required=False, default="Unknown Protein")
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
    parameters = ReferenceField(ParametersDBModel, required=True, reverse_delete_rule=CASCADE)
    log = EmbeddedDocumentListField(LogEntryDBModel, required=False)
    results = ReferenceField(ResultsDBModel, required=False, default=ResultsDBModel().save(), reverse_delete_rule=CASCADE)

    meta = {'collection': 'job', 'queryset_class': LoggerQuerySet}
