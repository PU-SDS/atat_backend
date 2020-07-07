import os
import subprocess
from io import StringIO
from os import path, rename
from typing import Tuple

from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from app.preprocess import PreProcessor
from flask import current_app as app


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

            SPLIT_MSA_FILES_PREFIX: The suffix used for file naming when the trimmed,
            and aligned MSA file is split into HOST and RESERVOIR files.

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
    SPLIT_MSA_FILES_PREFIX = 'SEQUENCE_'
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
        self.SPLIT_MSA_FILES_PREFIX = path.join(self.JOB_LOCAL_PATH,
                                                app.config['SPLIT_SEQS_SUBFOLDER'],
                                                self.SPLIT_MSA_FILES_PREFIX)

    def align(self):
        """
            Multiple Sequence Alignment method using MAFFT

            Carries out MSA alignment using MAFFT.
        """

        # Use BioPython cli wrapper for MAFFT
        mafft_cline = MafftCommandline(app.config['MAFFT'],
                                       input=path.join(self.JOB_LOCAL_PATH,
                                                       PreProcessor.PPROCESSED_OUTPUT_FILENAME))

        # Extract the stdout. This usually throw errors when sequence is malformed
        stdout = mafft_cline()[0]

        # Read the sequences as fasta and write to local path
        alignment = AlignIO.read(StringIO(stdout), 'fasta')
        AlignIO.write(alignment, path.join(self.JOB_LOCAL_PATH, self.ALIGNED_MSA_OUTPUT_FILENAME), 'fasta')


    def trim(self):
        """
            Alignment trimming using TRIMAL

            Trims the MSA using the user provided gap penalty (default = 0.05).
        """

        # Build the arguments for TRIMAL
        trimal_args = ['-keepheader', '-in', path.join(self.JOB_LOCAL_PATH, self.ALIGNED_MSA_OUTPUT_FILENAME),
                       '-out', path.join(self.JOB_LOCAL_PATH, self.TRIMMED_MSA_OUTPUT_FILENAME), '-gt',
                       f'{self.trimal_gt}']

        # Use Python subprocess to call TRIMAL
        trimal = self._run_process([app.config['TRIMAL']] + trimal_args)

        # Extract the stderr and the exitcode
        stderr = trimal[0]
        return_code = trimal[1]

        # If exitcode is not 0, throw an error
        if return_code:
            raise AlignmentUnknownException(f'An unknown error occurred during the sequence '
                                            f'trimming phase using TRIMAL: {stderr}')

    def split(self):
        """
            Alignment file splitting.

            Splits the alignment file into HOST and RESERVOIR files.
        """

        # Build the arguments for csplit
        csplit_args = ['-z', path.join(self.JOB_LOCAL_PATH, self.TRIMMED_MSA_OUTPUT_FILENAME),
                       '-f', self.SPLIT_MSA_FILES_PREFIX, '/HOST/']

        # Use Python subprocess to call csplit
        csplit = self._run_process(['csplit'] + csplit_args)
        stderr = csplit[0]
        return_code = csplit[1]

        # If exitcode is not 0, throw an error
        if return_code:
            raise AlignmentUnknownException(f'An unknown error occurred during the sequence '
                                            f'splitting phase using csplit: {stderr}')

        # Rename the split files
        self._rename_split_seq_files(self.SPLIT_MSA_FILES_PREFIX)

    @classmethod
    def _rename_split_seq_files(cls, split_file_prefix):
        """
            Rename the split file created by  csplit
        """

        host_seq = f'{split_file_prefix}00'
        reservoir_seq = f'{split_file_prefix}01'

        if not path.isfile(host_seq) or not path.isfile(reservoir_seq):
            raise AlignmentUnknownException('The split sequence files from csplit are not present.')

        rename(host_seq, f'{split_file_prefix}HOST.fasta')
        rename(reservoir_seq, f'{split_file_prefix}RESERVOIR.fasta')

    @classmethod
    def _run_process(cls, args) -> Tuple[bytes, int]:
        """
            Runs processes and returns the stderr and exitcode.

            Returns:
                Tuple[bytes, int]: A tuple containing the stderr and exitcode.
        """

        process = subprocess.Popen(args, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stderr = process.communicate()[1]

        return stderr, process.returncode
