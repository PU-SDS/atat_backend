from datetime import datetime
from typing import List

from pydantic import BaseModel, Field, EmailStr


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
