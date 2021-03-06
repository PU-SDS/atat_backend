from datetime import datetime
from typing import List

from pydantic import BaseModel, Field, EmailStr

from ...core.models.models import MotifClasses


class MotifTransmission(BaseModel):
    position: int = Field(..., title='The position at which the transmission occurred.')
    sequence: str = Field(..., title='The amino-acid sequence of the variant that got transmitted.')
    source: MotifClasses = Field(..., title='The classification of the variant in the host dataset.')
    target: MotifClasses = Field(..., title='The classification of the variant in the reservoir dataset.')
    source_incidence: float = Field(..., description="The incidence in the source.")
    target_incidence: float = Field(..., description="The incidence in the target.")


class MotifTransmissions(BaseModel):
    transmissions: List[MotifTransmission] = Field(..., title='List of motif transmissions that were observed.')


class Variant(BaseModel):
    class Config:
        allow_population_by_field_name = True

    sequence: str = Field(..., title='The amino-acid sequence of the variant.')
    count: int = Field(..., title='The number of times the variant was observed (frequency).')
    incidence: float = Field(..., title='The incidence of the variant.')
    motif_long: MotifClasses = Field(..., title='The class the variant belongs to.', alias='motif')
    metadata: dict = Field(..., title='Header metadata derived from the FASTA headers.')


class Position(BaseModel):
    position: int = Field(..., title='The kmer position.')
    support: int = Field(..., title='The number of supports that the kmer position has.')
    variants: List[Variant] = Field(None, title='Kmer variants seen at the kmer position.')


class GroupedPosition(BaseModel):
    """
    The API model for a grouped position response.
    """

    dataset_one: Position = Field(..., title='The dataset one kmer position.')
    dataset_two: Position = Field(..., title='The dataset two kmer position.')


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


class Parameters(BaseModel):
    kmer_length: int = Field(..., title='The length of the kmers that need to be generated using DiMA.')
    header_format: List[str] = Field(..., title='The ordered list of components of the FASTA header.')
    protein_name: str = Field('Unknown Protein', title='The name of the protein being worked with.')
    email: EmailStr = Field(None, title='The email address of the user.')


class CreateVivaJobRequest(BaseModel):
    id: str = Field(..., title='The DIM job id.')
    dataset_one_positions: List[Position] = Field(..., title='DiMA positions for the host dataset.')
    dataset_two_positions: List[Position] = Field(..., title='DiMA positions for the reservoir dataset.')
    parameters: Parameters = Field(..., title='The parameters to be used to run the analysis.')


class CreateStandaloneJobRequest(BaseModel):
    """
    This is the model for a job submit request.
    """

    dataset_one: str = Field(..., title='The aligned FASTA from the dataset one.')
    dataset_two: str = Field(..., title='The aligned FASTA from the dataset two.')
    parameters: Parameters = Field(..., title='The parameters to be used to run the analysis.')
