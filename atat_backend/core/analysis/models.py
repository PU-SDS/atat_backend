from pydantic import BaseModel, Field
from ..models.models import MotifClasses


class MotifTransmission(BaseModel):
    class Config:
        allow_population_by_field_name = True
        use_enum_values = True

    position: int = Field(..., description="The kmer position at which this switch takes place")
    sequence: str = Field(..., description="The amino-acid sequence of the kmer.")
    source_incidence: float = Field(..., description="The incidence in the source.")
    target_incidence: float = Field(..., description="The incidence in the target.")
    source: MotifClasses = Field(..., description="The original (source) motif name.")
    target: MotifClasses = Field(..., description="The target motif name.")
