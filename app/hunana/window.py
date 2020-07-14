from itertools import islice
from Bio import SeqIO


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

    def __init__(self, seqs, kmer_len=9):
        """
            Args:
                seqs: An iterable of type SeqRecord.
                kmer_len: The length of the sliding window (Default: 9).
        """

        self.seqs = seqs
        self.kmer_len = kmer_len

    def kmers(self):
        """
            Extracts k-mers from all the sequences provided.

            Returns:
                Zip: A Zip object containing all k-mers.
        """

        seqs = map(lambda s: str(s.seq), self.seqs)
        kmers_seqs = map(lambda s: self._window(s, self.kmer_len), seqs)
        return zip(*kmers_seqs)

    @classmethod
    def _window(cls, seq, kmer_len):
        """
            The actual logic of sliding window.
        """

        it = iter(seq)
        result = ''.join(islice(it, kmer_len))
        if len(result) == kmer_len:
            yield result
        for elem in it:
            result = result[1:] + elem
            yield result
