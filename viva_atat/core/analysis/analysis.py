from functools import cached_property
from io import StringIO
from itertools import product
from typing import List

from dima import Dima
from dima.helpers import Results, Variant, Position

from .models import MotifTransmission


class Analysis(object):
    def __init__(self, kmer_length: int, header_format: List[str], host_sequences: str, reservoir_sequences: str):
        """
        This class does all the analysis for the ATAT module. It is used by the Celery task.

        :param kmer_length: The length of kmers to generate via DiMA.
        :param header_format: The format of the header as provided by the user.
        :param host_sequences: The FASTA sequences for the host (ex: Human)
        :param reservoir_sequences: The FASTA sequences for the reservoir

        :type kmer_length: int
        :type header_format: List[str]
        :type host_sequences: str
        :type reservoir_sequences: str
        """

        self.kmer_length = kmer_length
        self.header_format = header_format
        self.host_sequences = host_sequences
        self.reservoir_sequences = reservoir_sequences

    @cached_property
    def dima_host(self) -> Results:
        """
        Runs DiMA for the host sequences and caches the results.
        """

        return Dima(
            StringIO(self.host_sequences),
            self.kmer_length,
            header_fillna="Unknown",
            header_format='|'.join(self.header_format),
        ).run()

    @cached_property
    def dima_reservoir(self) -> Results:
        """
        Runs DiMA for the reservoir sequences and caches the results.
        """

        return Dima(
            StringIO(self.reservoir_sequences),
            self.kmer_length,
            header_fillna="Unknown",
            header_format='|'.join(self.header_format),
        ).run()

    def atat(self) -> List[MotifTransmission]:
        """
        Runs transmission analysis where we look at positions where a kmer sequences goes from being classified as
        a particular class (ex: Index) in the reservoir to being classified differently in the zoonotic host
        """

        zipped_positions = zip(self.dima_host.results, self.dima_reservoir.results)
        switches = list()

        for host_position, reservoir_position in zipped_positions:  # type: Position
            host_variants = host_position.variants  # type: List[Variant]
            reservoir_variants = reservoir_position.variants  # type: List[Variant]

            for host_variant, reservoir_variant in product(host_variants, reservoir_variants):  # type: Variant
                if host_variant.sequence == reservoir_variant.sequence:
                    if host_variant.motif_short != reservoir_variant.motif_short:
                        switches.append(
                            MotifTransmission(
                                position=host_position.position,
                                sequence=host_variant.sequence,
                                fromx=reservoir_variant.motif_long,
                                to=host_variant.motif_long,
                            )
                        )

        return switches
