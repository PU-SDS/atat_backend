from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from io import StringIO
from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
import os

def sequence_removeDuplicates(sequences="", min_length=0, por_n=100):
        # Create our hash table to add the sequences
        cleanSequences=[]
        hashSequences={}
        # Using the Biopython fasta parse we can read our fasta input
        for seq_record in sequences:
                # Take the current sequence
                sequence = str(seq_record.seq).upper()
                # Check if the current sequence is according to the user parameters
                if (len(sequence) >= min_length and (float(sequence.count("N"))/float(len(sequence)))*100 <= por_n):
                # If the sequence passed in the test "is it clean?" and it isn't in the
                # hash table, the sequence and its id are going to be in the hash
                        if sequence not in hashSequences:
                                cleanSequences.append(seq_record)
                                hashSequences[sequence] = seq_record.id
        return cleanSequences

def sequence_alignment(file=''):
        mafft_cline = MafftCommandline('/home/centos/vada-single/bioapps/mafft-linux64/mafft.bat',input=file)
        stdout, stderr = mafft_cline()
        align = AlignIO.read(StringIO(stdout), "fasta")
        return align

def sequence_alignmentqc(infile='',outfile='',gap='0.05'):
        status = os.system('/home/centos/vada-single/bioapps/trimal -keepheader -in '+infile+' -out '+outfile+' -gt '+gap)
        return status

def sequence_split(infile='',prefix='',host='',baseFolder=''):
        os.system("csplit -z "+infile+" -f "+baseFolder+prefix+" '/tag=Host/'")
        os.rename(baseFolder+prefix+"00",baseFolder+prefix+"_Host.fasta")
        os.rename(baseFolder+prefix+"01",baseFolder+prefix+"_Reservoir.fasta")

def fasta2clustal(infile='',outfile=''):
	align = AlignIO.read(infile,'fasta')
	AlignIO.write(align,outfile,"clustal")


def preprocess(app,analysisId="atat",hostSequenceFile="",reservoirSequenceFile="",removeDuplicates=True,gapThreshold=""):
	baseFolder = app.config['UPLOADS_DEFAULT_DEST']
	hostInputSequenceFile = baseFolder + "hostInputSequence/" + hostSequenceFile
	reservoirInputSequenceFile = baseFolder + "reservoirInputSequence/" + reservoirSequenceFile

	# Tag the sequence

	hostSequence=[]

	with open(hostInputSequenceFile, "r") as handle:
		for record in SeqIO.parse(handle, "fasta") :
			record.id=record.id+"|tag=Host"
			hostSequence.append(record)

	reservoirSequence=[]
	with open(reservoirInputSequenceFile, "r") as handle:
        	for record in SeqIO.parse(handle, "fasta") :
            		record.id=record.id+"|tag=Reservoir"
            		reservoirSequence.append(record)

	# Merge the sequences
	listOfAllSequences=hostSequence+reservoirSequence

	# remove duplicates if necessary
	#cleanSequenceList=sequence_removeDuplicates(sequences=listOfAllSequences)
	#if not removeDuplicates:
	cleanSequenceList=listOfAllSequences

	# Aligned the file
	cleanSequenceFile=app.config['UPLOADS_DEFAULT_DEST']+analysisId+".fasta"
	alignedSequenceFile=app.config['UPLOADS_DEFAULT_DEST']+analysisId+"_aligned.fasta"
	SeqIO.write(cleanSequenceList,cleanSequenceFile,"fasta")
	aligned_sequence=sequence_alignment(file=cleanSequenceFile)
	AlignIO.write(aligned_sequence,alignedSequenceFile,"fasta")

	# Alignment QC
	qcSequenceFile=app.config['UPLOADS_DEFAULT_DEST']+analysisId+"_qced.fasta"
	sequence_alignmentqc(infile=alignedSequenceFile,outfile=qcSequenceFile)

	# Split Sequence
	sequence_split(baseFolder=app.config['UPLOADS_DEFAULT_DEST'],infile=qcSequenceFile,prefix=analysisId)
	return app.config['UPLOADS_DEFAULT_DEST']+analysisId+"_Host.fasta",app.config['UPLOADS_DEFAULT_DEST']+analysisId+"_Reservoir.fasta"
    #fasta2clustal(infile=app.config['UPLOADS_DEFAULT_DEST']+analysisId+"_Animal.fasta",outfile=app.config['UPLOADS_DEFAULT_DEST']+analysisId+"_Animal.aln")
	#fasta2clustal(infile=app.config['UPLOADS_DEFAULT_DEST']+analysisId+"_Input.fasta",outfile=app.config['UPLOADS_DEFAULT_DEST']+analysisId+"_Input.aln")
