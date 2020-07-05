import subprocess
from io import StringIO
from os import path
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
            align(object)
            trim(object)

        Constants:
            ALIGNED_MSA_OUTPUT_FILENAME: The filename used for the MAFFT aligned
            sequence file.

            TRIMMED_MSA_OUTPUT_FILENAME: The filename used for the TRIMAI trimmed
            sequence file.

            JOB_LOCAL_PATH: The local path to the job. This is built automatically
            when the class is instantiated using the job id and the common job path
            defined in the configuration file.

        Example:
            >>> from app.preprocess import Alignment
            >>> new_alignment = Alignment('id_5165432', 0.05)
            >>> new_alignment.align()
            >>> new_alignment.trim()
    """

    ALIGNED_MSA_OUTPUT_FILENAME = 'aligned.fasta'
    TRIMMED_MSA_OUTPUT_FILENAME = 'aligned_trimmed.fasta'
    JOB_LOCAL_PATH = ''

    def __init__(self, job_id: str, trimal_gt: float = 0.05):
        """
            Args:
                job_id (str): The current job ID.
                trimal_gt (float): The gap threshold for TRIMAL.
        """
        self.job_id = job_id
        self.trimal_gt = trimal_gt
        self.JOB_LOCAL_PATH = path.join(app.config['JOBS_FOLDER'], self.job_id)

    def align(self):
        """
            Multiple Sequence Alignment method using MAFFT

            Carries out MSA alignment using MAFFT.
        """

        mafft_cline = MafftCommandline(app.config['MAFFT'],
                                       input=path.join(self.JOB_LOCAL_PATH, PreProcessor.OUTPUT_FILENAME))

        try:
            std_out, std_err = mafft_cline()
        except ApplicationError as error:
            raise AlignmentUnknownException(f'An unknown error occurred during the sequence '
                                            f'alignment phase using MAFFT: {std_err}')

        alignment = AlignIO.read(StringIO(std_out), 'fasta')
        AlignIO.write(alignment, path.join(self.JOB_LOCAL_PATH, self.ALIGNED_MSA_OUTPUT_FILENAME))

    def trim(self):
        """
            Alignment trimming method using TRIMAL

            Trims the MSA using the user provided gap penalty (default = 0.05).
        """
        trimal_args = ['-in', path.join(self.jobpath, self.ALIGNED_MSA_OUTPUT_FILENAME),
                       '-out', path.join(self.jobpath, self.TRIMMED_MSA_OUTPUT_FILENAME), '-gt', self.trimal_gt]

        trimal = subprocess.Popen([app.config['TRIMAL']] + trimal_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        std_out, std_err = trimal.communicate()
        return_code = trimal.returncode

        if return_code:
            raise AlignmentUnknownException(f'An unknown error occurred during the sequence '
                                            f'trimming phase using TRIMAL: {std_err}')
