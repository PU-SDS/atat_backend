import itertools
from os import path
from typing import Generator, Iterable
from Bio import SeqIO, SeqRecord
from Bio.Alphabet import generic_protein
from flask import app


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

            JOB_LOCAL_PATH: The local path to the job. This is built automatically
            when the class is instantiated using the job id and the common job path
            defined in the configuration file.

        Example:
            >>> from app.preprocess import PreProcessor
            >>> preprocessor = PreProcessor('id_5165432')
            >>> preprocessor.preprocess()
    """

    INPUT_HOST_FILENAME = 'host.fasta'
    INPUT_RESERVOIR_FILENAME = 'reservoir.fasta'
    PPROCESSED_OUTPUT_FILENAME = 'processed.fasta'
    JOB_LOCAL_PATH = ''

    def __init__(self, job_id: str):
        """
        Args:
             job_id (str): The current job ID.
        """
        self.job_id = job_id
        self.JOB_LOCAL_PATH = path.join(app.config['JOBS_FOLDER'], self.job_id)

    def preprocess(self):
        """
            Tags the sequences accordingly (HOST, RESERVOIR), merges the sequences
            into a single file and removes any duplicate sequences.
        """

        # Tags the sequences based on origin (HOST, RESERVOIR)
        reservoir = path.join(self.JOB_LOCAL_PATH, self.INPUT_RESERVOIR_FILENAME)
        host = path.join(self.JOB_LOCAL_PATH, self.INPUT_HOST_FILENAME)

        host_tagged = self._read_tag_seq(host, 'HOST')
        reservoir_tagged = self._read_tag_seq(reservoir, 'RESERVOIR')

        # Merge the host and reservoir sequences
        merged_seqs = itertools.chain(reservoir_tagged, host_tagged)

        if not is_iterable_of_type(merged_seqs, SeqRecord):
            raise PreprocessorException('Merging of sequences failed.')

        # Gets the unique sequences
        unique_seq_records = self._get_unique_sequences(merged_seqs)

        # Writes merged entries into a file
        SeqIO.write(unique_seq_records, path.join(self.JOB_LOCAL_PATH, self.PPROCESSED_OUTPUT_FILENAME))

    @classmethod
    def _read_tag_seq(cls, seq_path, tag) -> Iterable[SeqRecord]:
        """
            Reads a file of sequences and tags
        """

        with open(seq_path, 'r') as f:
            for record in SeqIO.parse(f, 'fasta', generic_protein):
                record.id += f'|{tag}'
                yield record

    @classmethod
    def _get_unique_sequences(cls, sequences) -> Iterable[SeqRecord]:
        """
            Gets a list of unique sequences.
        """

        hashed_seqs = set()

        for seq_record in sequences:
            seq = seq_record.seq

            if seq not in hashed_seqs:
                hashed_seqs.add(seq)
                yield seq_record


def is_iterable_of_type(itera, klass):
    return all(map(lambda x: isinstance(x, klass), itera))
    # for x in itera:
    #     if not isinstance(x, klass):
    #         return False
    # return True
