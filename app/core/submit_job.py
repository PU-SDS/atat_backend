from os.path import join
from flask import current_app as app
from Bio import SeqIO
from Bio.Alphabet import IUPAC

from app.hunana import SlidingWindow
from app.hunana.entropy import NormalizedEntropy
from app.core import Constants
from app.warehousing.mongodb.write import MongoDBWrite
from app.warehousing.mongodb.constants import JobStatuses
from app.metadata  import UpdateMetadata
from app.notifications import EmailNotification


class SubmitJob(object):
    def __init__(self, job_id, email, kmer_length=9, alignment=False):
        """
            Submits the job to the sub queue
            :param job_id: The job id
            :param alignment: Whether we need MSA for the input sequences

            :type job_id: str
            :type alignment: bool
        """

        self.job_id = job_id
        self.email= email
        self.kmer_length = kmer_length
        self.alignment = alignment

    def run(self):
        # Build the job folder path
        job_folder = join(app.config['JOBS_FOLDER'], self.job_id)

        # Warehousing instance
        insert_id = db_write = MongoDBWrite(self.job_id)

        # Update job status to PROCESSING
        db_write.update_job_status(JobStatuses.JOB_PROCESSING)

        # Read the host and reservoir sequence files
        host_sequence_file = open(join(job_folder, Constants.HOST_SEQUENCE_FILENAME))
        reservoir_sequence_file = open(join(job_folder, Constants.RESERVOIR_SEQUENCE_FILENAME))

        # Parse the host and reservoir sequence files
        host_sequence_parsed = SeqIO.parse(host_sequence_file, 'fasta', IUPAC.protein)
        reservoir_sequence_parsed = SeqIO.parse(reservoir_sequence_file, 'fasta', IUPAC.protein)

        # Get the correct metadata for the headers
        host_updated_metadata = UpdateMetadata(sequences=host_sequence_parsed).update()
        reservoir_updated_metadata = UpdateMetadata(sequences=reservoir_sequence_parsed).update()

        # Get k-mer sequences from both host and reservoir sequences
        host_kmers = list(SlidingWindow(host_updated_metadata, self.kmer_length).run())
        reservoir_kmers = list(SlidingWindow(reservoir_updated_metadata, self.kmer_length).run())

        # Calculate the entropy for each k-mer position of both host and sequence
        host_kmer_positions = NormalizedEntropy(10000, 10, host_kmers).run()
        reservoir_kmer_positions = NormalizedEntropy(10000, 10, reservoir_kmers).run()

        # Data warehousing
        host_data = list(host_kmer_positions)
        reservoir_data = list(reservoir_kmer_positions)

        # Write results to the database
        db_write_status = db_write.save_positions(host_data, reservoir_data)

        # Close the file handles
        host_sequence_file.close()
        reservoir_sequence_file.close()

        # Update job status to COMPLETED
        db_write.update_job_status(JobStatuses.JOB_COMPLETED)

        # Send email notification
        EmailNotification(self.email, self.job_id).send()

        # Print the status of the database write
        print(db_write_status)
