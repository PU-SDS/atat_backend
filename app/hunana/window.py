from collections import Counter
from itertools import islice

from app.hunana.position import Position, Variant


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

        self.seqs = seqs
        self.kmer_len = kmer_len

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

    @classmethod
    def _variant_counter(cls, kmers):
        """
            Counts the number of various variants at each kmer position

            Returns:
                A Generator object for a counter for all variants at each kmer position
        """

        for nom in kmers:
            counter = Counter(nom)
            del counter['illegal-char']
            yield counter

    def _create_position_objects(self, counters):
        """
            Create kmer position objects with entropy, position and variant properties

            Returns:
                A Generator object for Position objects for each kmer position
        """

        for position, counter in enumerate(counters):
            position_object = Position()
            position_object.POSITION = position
            position_object.VARIANTS = self._create_variant_objects(counter)
            position_object.VARIANTS_FLATTENED = counter.elements()

            yield position_object

    def _create_variant_objects(self, variant_counter):
        variant_count = sum(variant_counter.values())
        variant_objects = map(lambda x: Variant(x[0], x[1], self._calc_conservation(x[1], variant_count)),
                              variant_counter.items())
        return list(variant_objects)

    @classmethod
    def _calc_conservation(cls, variant_hits, total_hits):
        conservation = (float(variant_hits) * 100) / float(total_hits)
        return conservation

    def run(self):
        kmers = self._kmers()
        variant_counters = self._variant_counter(kmers)
        kmer_position_objects = self._create_position_objects(variant_counters)
        return kmer_position_objects
