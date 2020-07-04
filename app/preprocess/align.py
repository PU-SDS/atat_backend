import subprocess
from io import StringIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MafftCommandline
from Bio.Application import ApplicationError
from flask import app


class AlignmentUnknownException(Exception):
    def __int__(self, error):
        self.error = error
        super(AlignmentUnknownException, self).__init__(self.error)


class Alignment(object):

    """
    Alignment Class

    Contains all logic related to MSA including trimming and alignment itself.

    Functions:
        align(object, file: str) -> str
        trim(object, file: str) -> str

    Example:
        new_alignment = Alignment(file)
        new_alignment.align()
        new_alignment.trim()
    """

    def __init__(self, basefile: str, msa=None):
        self.basefile = basefile
        self.msa = msa

    def align(self) -> str:
        """
            Multiple Sequence Alignment Method using MAFFT

            This method carries out MSA alignment using MAFFT.

            Returns:
                str: Returns a string containing a MSA
        """
        mafft_cline = MafftCommandline(app.config['MAFFT'], input=self.basefile)

        try:
            stdout, stderr = mafft_cline()
            self.msa = stdout
        except ApplicationError as error:
            raise AlignmentUnknownException(f'An unknown error occurred during the sequence '
                                            f'alignment phase using MAFFT: {error}')
        self.msa = AlignIO.read(StringIO(self.msa), 'fasta')

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
