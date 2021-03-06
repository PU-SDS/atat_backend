from io import StringIO
from typing import List

from dima import Dima
from dima.helpers import Results

from .models import MotifTransmission


class Analyses(object):
    @staticmethod
    def dima_analysis(sequences: str, kmer_length: int, header_format: List[str]) -> Results:
        """
        A common method to run all DiMA analyses within the context of ATAT.

        :param sequences: The FASTA sequences for the host (ex: Human)
        :param kmer_length: The length of kmers to generate via DiMA.
        :param header_format: The format of the header as provided by the user.

        :type sequences: str
        :type kmer_length: int
        :type header_format: List[str]

        :returns: A Results object that contains the full DiMA analysis results.
        """

        return Dima(
            StringIO(sequences),
            kmer_length,
            header_fillna="Unknown",
            header_format='|'.join(header_format),
        ).run()

    @staticmethod
    def at_analysis(host_dima_positions: List[dict], reservoir_dima_positions: List[dict]) -> List[MotifTransmission]:
        """
        Runs transmission analysis where we look at positions where a kmer sequences goes from being classified as
        a particular class (ex: Index) in the reservoir to being classified differently in the zoonotic host.

        :param host_dima_positions: Kmer positions generated by DiMA for the source sequences.
        :param reservoir_dima_positions: Kmer positions generated by DiMA for the reservoir sequences.

        :type host_dima_positions: List[dict]
        :type reservoir_dima_positions: List[dict]

        :returns: A list of MotifTransmission objects for each motif transmission observed.
        """

        zipped_positions = zip(host_dima_positions, reservoir_dima_positions)
        switches = list()

        for host_position, reservoir_position in zipped_positions:  # type: dict
            host_variants = host_position.get('variants')
            reservoir_variants = reservoir_position.get('variants')

            if not host_variants or not reservoir_variants:
                continue

            for host_variant in host_variants:
                reservoir_variant_match = [reservoir_variant for reservoir_variant in reservoir_variants if
                                           reservoir_variant.get('sequence') == host_variant.get('sequence')]

                if not reservoir_variant_match:
                    continue

                if host_variant.get('motif_long') != reservoir_variant_match[0].get('motif_long'):
                    switches.append(
                        MotifTransmission(
                            position=host_position.get('position'),
                            sequence=host_variant.get('sequence'),
                            source=host_variant.get('motif_long'),
                            target=reservoir_variant_match[0].get('motif_long'),
                            source_incidence=host_variant.get('incidence'),
                            target_incidence=reservoir_variant_match[0].get('incidence')
                        )
                    )

        return switches
