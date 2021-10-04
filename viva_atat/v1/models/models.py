from datetime import datetime
from typing import List

from pydantic import BaseModel, Field, EmailStr

from ...core.models.models import MotifClasses


class MotifTransmission(BaseModel):
    position: int = Field(..., title='The position at which the transmission occurred.')
    sequence: str = Field(..., title='The amino-acid sequence of the variant that got transmitted.')
    source: MotifClasses = Field(..., title='The classification of the variant in the host dataset.')
    target: MotifClasses = Field(..., title='The classification of the variant in the reservoir dataset.')


class MotifTransmissions(BaseModel):
    transmissions: List[MotifTransmission] = Field(..., title='List of motif transmissions that were observed.')


class Parameters(BaseModel):
    kmer_length: int = Field(..., title='The length of the kmers generated by DiMA.')
    header_format: List[str] = Field(..., title='The ordered list of components of the FASTA headers.')


class Variant(BaseModel):
    sequence: str = Field(..., title='The amino-acid sequence of the variant.')
    count: int = Field(..., title='The number of times the variant was observed (frequency).')
    incidence: float = Field(..., title='The incidence of the variant.')
    motif: MotifClasses = Field(..., title='The class the variant belongs to.')
    metadata: dict = Field(..., title='Header metadata derived from the FASTA headers.')


class Position(BaseModel):
    position: int = Field(..., title='The kmer position.')
    supports: int = Field(..., title='The number of supports that the kmer position has.')
    variants: List[Variant] = Field(..., title='Kmer variants seen at the kmer position.')


class GroupedPosition(BaseModel):
    """
    The API model for a grouped position response.
    """

    host: Position = Field(..., title='The host kmer position.')
    reservoir: Position = Field(..., title='The reservoir kmer position.')


class JobLogEntry(BaseModel):
    """
    The API model for a single log entry for a job.
    """

    id: str = Field(..., title='The id of the log entry.')
    flag: str = Field(..., title='The flag of the log entry indicating severity.')
    timestamp: datetime = Field(..., title='The time at which the log entry was produced.')
    message: str = Field(..., title='The message body of the log entry.')


class JobLogs(BaseModel):
    """
    The API model for a list of log entries of a job.
    """

    logs: List[JobLogEntry] = Field(..., title='List of log entries for the job.')


class CreateJobParameters(BaseModel):
    kmer_length: int = Field(..., title='The length of the kmers that need to be generated using DiMA.')
    header_format: List[str] = Field(..., title='The ordered list of components of the FASTA header.')
    email: EmailStr = Field(None, title='The email address of the user.')


class CreateJobRequest(BaseModel):
    """
    This is the model for a job submit request.
    """

    host_sequences: str = Field(..., title='The aligned FASTA from the host species.')
    reservoir_sequences: str = Field(..., title='The aligned FASTA from the reservoir species.')
    parameters: CreateJobParameters = Field(..., title='The parameters to be used to run the analysis.')
