import itertools
from os import path
from typing import Generator, Iterable

from Bio import SeqIO, SeqRecord
from Bio.Alphabet import generic_protein


class PreprocessorException(Exception):
    def __init__(self, msg):
        self.msg = msg
        super(PreProcessor, self).__init__(self.msg)


class PreProcessor(object):
    """
        Pre-processing Class

        Performs all critical pre-processing such as appropriate tagging of sequences,
        merging and removal of duplicates.

        Constants:
            INPUT_HOST_FILENAME: The filename of the host sequence file.

            INPUT_RESERVOIR_FILENAME: The filename of the reservoir sequence file.

            OUTPUT_FILENAME: The filenamed used for the fully pre-processed sequence file.

        Example:
            >>> from app.preprocess import PreProcessor
            >>> preprocessor = PreProcessor('/jobs/id_5165432')
            >>> preprocessor.preprocess()
    """

    INPUT_HOST_FILENAME = 'host.fasta'
    INPUT_RESERVOIR_FILENAME = 'reservoir.fasta'
    OUTPUT_FILENAME = 'processed.fasta'

    def __init__(self, jobpath: str):
        """
        Args:
             jobpath (str): The folder path to the job.
        """
        self.jobpath = jobpath

    def preprocess(self):
        """
            Tags the sequences accordingly (HOST, RESERVOIR), merges the sequences
            into a single file and removes any duplicate sequences.
        """

        # Tags the sequences based on origin (HOST, RESERVOIR)
        reservoir = path.join(self.jobpath, self.INPUT_RESERVOIR_FILENAME)
        host = path.join(self.jobpath, self.INPUT_HOST_FILENAME)

        reservoir_tagged = self._read_seq(host, 'HOST')
        host_tagged = self._read_seq(reservoir, 'RESERVOIR')

        # Merges the host and reservoir sequences into a single file
        if not is_iterable_of_type(itertools.chain(reservoir_tagged, host_tagged), SeqRecord):
            raise PreprocessorException('Merging of sequences failed.')
        unique_seq_records = self._get_unique_sequence(reservoir_tagged + host_tagged)

        SeqIO.write(unique_seq_records, path.join(self.jobpath, self.OUTPUT_FILENAME))

    @classmethod
    def _read_seq(cls, path, tag) -> Iterable[SeqRecord]:
        """
            Reads a file of sequences and tags
        """

        with open(path, 'r') as f:
            for record in SeqIO.parse(f, 'fasta', generic_protein):
                record.id += f'|{tag}'
                yield record

    @classmethod
    def _get_unique_sequence(cls, sequences) -> Iterable[SeqRecord]:
        """
            Gets a list of unique sequences.
        """

        hashed_seqs = set()

        for seq_record in sequences:
            seq = seq_record.seq

            if seq not in hashed_seqs:
                hashed_seqs.add(seq)
                yield seq_record


def is_iterable_of_type(iter, klass):
    # return all(map(lambda x: isinstance(x, klass), iter))
    for x in iter:
        if not isinstance(x, klass):
            return False
    return True
