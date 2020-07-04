import subprocess
from io import StringIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MafftCommandline
from Bio.Application import ApplicationError
from flask import app
from app.preprocess import PreProcessor


class AlignmentUnknownException(Exception):
    def __int__(self, error):
        self.error = error
        super(AlignmentUnknownException, self).__init__(self.error)


class Alignment(object):
    """
        Alignment Class

        Contains all the logic related to MSA including trimming, splitting and alignment itself.

        Functions:
            align(object, jobpath: str)
            trim(object, file: str) -> str

        Constants:
            ALIGNED_MSA_OUTPUT_FILENAME: The filename used for the MAFFT aligned
            sequence file.

        Example:
            >>> from app.preprocess import Alignment
            >>> new_alignment = Alignment(jobpath)
            >>> new_alignment.align()
            >>> new_alignment.trim()
    """

    ALIGNED_MSA_OUTPUT_FILENAME = 'aligned.fasta'

    def __init__(self, jobpath: str, msa=None):
        self.jobpath = jobpath
        self.msa = msa

    def align(self) -> str:
        """
            Multiple Sequence Alignment Method using MAFFT

            This method carries out MSA alignment using MAFFT.
        """
        mafft_cline = MafftCommandline(app.config['MAFFT'],
                                       input=self.jobpath + f'/{PreProcessor.OUTPUT_FILENAME}')

        try:
            stdout, stderr = mafft_cline()
            self.msa = stdout
        except ApplicationError as error:
            raise AlignmentUnknownException(f'An unknown error occurred during the sequence '
                                            f'alignment phase using MAFFT: {error}')

        alignment = AlignIO.read(StringIO(self.msa), 'fasta')
        AlignIO.write(alignment, self.jobpath + f'/{self.ALIGNED_MSA_OUTPUT_FILENAME}')

    def trim(self, file: str):
        """
            Alignment trimming method using TRIMAL

            This method trims the MSA using the user provided gap penalty (default = 0.05).

            Parameters:
                file (str): The path to the aligned sequences file.

            Returns:
                A file path to the aligned and trimmed FASTA file.

        """
        if isinstance(self.msa, MultipleSeqAlignment):
            trimal = subprocess.Popen(
                [app.config['TRIMAL'], '']
            )
