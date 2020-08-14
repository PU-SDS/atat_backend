from collections import Counter
from itertools import islice

from app.hunana.position import Position, Variant, VariantDict


class SlidingWindow(object):
    """
        Sliding Window Class

        Gets a zip object containing k-mers from from all the sequences provided classified by position.

        Example:
            >>> from app.hunana import SlidingWindow
            >>> window = SlidingWindow(iter[SeqRecord], 9)
            >>> kmers = window.kmers()
            >>> kmer1 = list(kmers[0]) # First k-mer of all sequences
            >>> kmer2 = list(kmers[1]) # Second k-mer of all sequences
    """

    DISALLOWED_CHARS = {'-', 'X', 'B', 'J', 'Z', 'O', 'U'}

    def __init__(self, seqs, kmer_len=9):
        """
            Args:
                seqs: An iterable of type SeqRecord.
                kmer_len: The length of the sliding window (Default: 9).
        """

        self.seqs = list(seqs)
        self.kmer_len = kmer_len
        self.idx = [x.description for x in self.seqs]

    def _create_variant_dict(self, iterable):
        positions = []

        for c in iterable:
            unique = VariantDict(list)
            for a, b in enumerate(c):
                unique[b].append(a)
            unique.pop('illegal-char', None)
            positions.append(unique)
        return positions

    def _kmers(self):
        """
            Extracts k-mers from all the sequences provided.

            Returns:
                Zip: A Zip object containing all k-mers.
        """

        seqs = map(lambda s: str(s.seq), self.seqs)
        kmers_seqs = map(lambda s: self._window(s, self.kmer_len), seqs)
        return zip(*kmers_seqs)

    def _window(self, seq, kmer_len):
        """
            The actual logic of sliding window.

            Returns:
                A Generator object for all kmers within the alignment
        """

        it = iter(seq)
        result = ''.join(islice(it, kmer_len))
        if len(result) == kmer_len:
            if all(ele not in result for ele in self.DISALLOWED_CHARS):
                yield result
            else:
                yield 'illegal-char'

        for elem in it:
            result = result[1:] + elem
            if all(ele not in result for ele in self.DISALLOWED_CHARS):
                yield result
            else:
                yield 'illegal-char'

    def _create_position_objects(self, counters, variant_dict):
        """
            Create kmer position objects with entropy, position and variant properties

            Returns:
                A Generator object for Position objects for each kmer position
        """

        for idx, counter in enumerate(counters):
            yield Position(
                position=idx,
                sequences=self._create_variant_objects(idx, counter),
                variants_flattened=list(counter.elements()),
                variant_dict=variant_dict[idx]
            )

    def _create_variant_objects(self, idx, variant_counter: Counter) -> list:
        num_variants = sum(variant_counter.values())
        variant_objects = [
            Variant(idx, seq, count, self._calc_conservation(count, num_variants))
            for seq, count
            in variant_counter.items()
        ]
        return variant_objects

    @classmethod
    def _calc_conservation(cls, variant_hits, total_hits):
        conservation = (float(variant_hits) * 100) / float(total_hits)
        return float(conservation)

    def run(self):
        kmers = self._kmers()
        variant_dicts = self._create_variant_dict(kmers)
        variant_counters = [x.get_counter() for x in variant_dicts]
        kmer_position_objects = self._create_position_objects(variant_counters, variant_dicts)
        return kmer_position_objects
