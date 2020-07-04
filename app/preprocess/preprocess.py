from io import StringIO
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

        Parameters:
            jobpath (str): The folder path to the job.

        Functions:
            preprocess(object)

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
        self.jobpath = jobpath

    def preprocess(self):
        """
            Preprocessor method.

            This method tags the sequences accordingly (HOST, RESERVOIR), merges the sequences
            into a single file and removes any duplicate sequences.
        """

        # Tags the sequences based on origin (HOST, RESERVOIR)
        reservoir = self.jobpath + f'/{self.INPUT_RESERVOIR_FILENAME}'
        host = self.jobpath + f'/{self.INPUT_HOST_FILENAME}'

        reservoir_tagged = []
        host_tagged = []

        with open(reservoir, 'r') as reservoir_file, open(host, 'r') as host_file:
            for seq_record in SeqIO.parse(host_file, 'fasta', generic_protein):
                seq_record.id += '|HOST'
                host_tagged.append(seq_record)

            for seq_record in SeqIO.parse(reservoir_file, 'fasta', generic_protein):
                seq_record.id += '|RESERVOIR'
                reservoir_tagged.append(seq_record)

        # Merges the host and reservoir sequences into a single file
        if isinstance(reservoir_tagged, list(SeqRecord) and
                                        isinstance(host_tagged, list(SeqRecord))):
            merged_seqs: list[SeqRecord] = reservoir_tagged + host_tagged
        else:
            raise PreprocessorException('Merging of sequences failed.')

        # Removes any duplicate sequences
        hashed_seqs = {}
        unique_seq_records = []

        for seq_record in merged_seqs:
            seq = seq_record.seq

            if seq not in hashed_seqs:
                unique_seq_records.append(seq_record)
                hashed_seqs[seq] = seq_record.id

        SeqIO.write(unique_seq_records, self.jobpath + self.OUTPUT_FILENAME)
