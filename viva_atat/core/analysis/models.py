from pydantic import BaseModel, Field
from ..models.models import MotifClasses


class MotifTransmission(BaseModel):
    position: int = Field(..., description="The kmer position at which this switch takes place")
    sequence: str = Field(..., description="The amino-acid sequence of the kmer.")
    fromx: MotifClasses = Field(..., alias="from", description="The original (source) motif name.")
    to: MotifClasses = Field(..., description="The target motif name.")
