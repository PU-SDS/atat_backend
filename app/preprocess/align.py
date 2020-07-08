import subprocess
from io import StringIO
from os import path, rename
from typing import Tuple

from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from flask import current_app as app

from app.preprocess import PreProcessor


class AlignmentUnknownException(Exception):
    def __int__(self, error):
        self.error = error
        super(AlignmentUnknownException, self).__init__(self.error)


class Alignment(object):
    """
        Alignment Class

        Contains all the logic related to MSA including trimming, splitting, fasta2clustal and alignment itself.

        Functions:
            align(object)
            trim(object)
            split(object)
            toclustal(object)

        Constants:
            ALIGNED_MSA_OUTPUT_FILENAME: The filename used for the MAFFT aligned
            sequence file.

            TRIMMED_MSA_OUTPUT_FILENAME: The filename used for the TRIMAI trimmed
            sequence file.

            SPLIT_MSA_FILES_PREFIX: The suffix used for file naming when the trimmed,
            and aligned MSA file is split into HOST and RESERVOIR files. This contains
            the full path to the split files as the prefix.

            JOB_LOCAL_PATH: The local path to the job. This is built automatically
            when the class is instantiated using the job id and the common job path
            defined in the configuration file.

        Example:
            >>> from app.preprocess import Alignment
            >>> new_alignment = Alignment('id_5165432', 0.05)
            >>> new_alignment.align()
            >>> new_alignment.trim()
            >>> new_alignment.split()
            >>> new_alignment.toclustal()
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
        """

        processed_seqs_file = path.join(self.JOB_LOCAL_PATH,
                                        PreProcessor.PPROCESSED_OUTPUT_FILENAME)

        if not path.isfile(processed_seqs_file):
            raise AlignmentUnknownException('The processed sequence file is not present. '
                                            'Run the preprocessor first.')

        # Use BioPython cli wrapper for MAFFT
        mafft_cline = MafftCommandline(app.config['MAFFT'],
                                       input=path.join(self.JOB_LOCAL_PATH,
                                                       PreProcessor.PPROCESSED_OUTPUT_FILENAME))

        # Extract the stdout. This usually throw errors when sequence is malformed
        stdout = mafft_cline()[0]

        # Read the sequences as fasta and write to local path
        alignment = AlignIO.read(StringIO(stdout), 'fasta')
        AlignIO.write(alignment, path.join(self.JOB_LOCAL_PATH, self.ALIGNED_MSA_OUTPUT_FILENAME),
                      'fasta')

    def trim(self):
        """
            Trims the MSA using TRIMAL.
        """

        aligned_seqs_file = path.join(self.JOB_LOCAL_PATH, self.ALIGNED_MSA_OUTPUT_FILENAME)

        if not path.isfile(aligned_seqs_file):
            raise AlignmentUnknownException('The aligned sequence file is not present. '
                                            'Run align() first.')

        # Build the arguments for TRIMAL
        trimal_args = ['-keepheader', '-in', aligned_seqs_file,
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
            Splits the alignment file into HOST and RESERVOIR files.
        """

        trimmed_seqs_file = path.join(self.JOB_LOCAL_PATH, self.TRIMMED_MSA_OUTPUT_FILENAME)

        if not path.isfile(trimmed_seqs_file):
            raise AlignmentUnknownException('The trimmed sequence file is not present. Run trim() first.')

        # Build the arguments for csplit
        csplit_args = ['-z', trimmed_seqs_file,
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

    def toclustal(self):
        """
            Converts the split FASTA files into CLUSTAL format.
        """

        host_seqs_file = f'{self.SPLIT_MSA_FILES_PREFIX}HOST.fasta'
        reservoir_seqs_file = f'{self.SPLIT_MSA_FILES_PREFIX}RESERVOIR.fasta'

        if not path.isfile(host_seqs_file) or not path.isfile(reservoir_seqs_file):
            raise AlignmentUnknownException('The split sequences are not present. Run split() first.')

        # Read the host and reservoir FASTA files
        host_fasta_alignment = AlignIO.read(host_seqs_file, 'fasta')
        reservoir_fasta_alignment = AlignIO.read(reservoir_seqs_file, 'fasta')

        # Write the host and reservoir sequences as CLUSTAL format
        AlignIO.write(host_fasta_alignment, f'{self.SPLIT_MSA_FILES_PREFIX}HOST_CLUSTAL.txt')
        AlignIO.write(reservoir_fasta_alignment, f'{self.SPLIT_MSA_FILES_PREFIX}RESERVOIR_CLUSTAL.txt')

    @classmethod
    def _rename_split_seq_files(cls, split_file_prefix: str):
        """
            Renames the split files created by csplit
        """

        reservoir_seq = f'{split_file_prefix}00'
        host_seq = f'{split_file_prefix}01'

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
