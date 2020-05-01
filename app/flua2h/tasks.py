from flua2h.celery import app
from flua2h.atat import Atat
from common.sequence import sequence_split,sequence_alignmentqc,sequence_download,getTemporaryFileName,sequence_removeDuplicates,sequence_alignment
from common.jobmanagement import jobUpdate
from Bio import SeqIO
from Bio import AlignIO
from common.hunana import Hunana
import os
import json
import logging

def fileCleanUp(listOfFiles):
	for file in listOfFiles:
		os.remove(file)

@app.task(trail=True)
def flua2h(id="",index="",protein="",host=""):
	jobUpdate(id,"EXECUTING","Sequence Download")
	listOfAllSequences=[]
	animalSequence = sequence_download(index,protein,host,'Animal')
	humanSequence = sequence_download(index,protein,'Human','Human')
	listOfAllSequences=animalSequence+humanSequence

	baseFolder = "temp/"
	#baseFile = getTemporaryFileName(ext='')
	baseFile = id
	sequenceFile=baseFolder+baseFile+".fasta"
	alignedSequenceFile=baseFolder+"aligned_"+baseFile+".fasta"
	cleanSequenceFile=baseFolder+"clean_"+baseFile+".fasta"
	qcSequenceFile=baseFolder+"qc_"+baseFile+".fasta"

	jobUpdate(id,"EXECUTING","Sequence Alignment")
	cleanSequenceList=sequence_removeDuplicates(sequences=listOfAllSequences)
	SeqIO.write(cleanSequenceList,cleanSequenceFile,"fasta")
	aligned_sequence=sequence_alignment(file=cleanSequenceFile)
	AlignIO.write(aligned_sequence,alignedSequenceFile,"fasta")
	sequence_alignmentqc(infile=alignedSequenceFile,outfile=qcSequenceFile)
	sequence_split(baseFolder=baseFolder,infile=qcSequenceFile,prefix=baseFile,host=host)
	animalFastaFile=baseFolder+baseFile+"_Animal.fasta"
	humanFastaFile=baseFolder+baseFile+"_Human.fasta"

	listOfFiles=(alignedSequenceFile,cleanSequenceFile,qcSequenceFile,animalFastaFile,humanFastaFile)
	jobUpdate(id,"EXECUTING","HUNANA Analysis")
	try:
		HUNANA=Hunana()
		humanHunanaResult=HUNANA.run(id=baseFile,inputfile=humanFastaFile,protein=protein,source="human")
		animalHunanaResult=HUNANA.run(id=baseFile,inputfile=animalFastaFile,protein=protein,source="animal")
		outputstream=open(baseFolder+baseFile+"_animal.hunana","w")
		json.dump(animalHunanaResult,outputstream)
		outputstream.close()
		outputstream=open(baseFolder+baseFile+"_human.hunana","w")
		json.dump(humanHunanaResult,outputstream)
		outputstream.close()
	except Exception as error:
		jobUpdate(id,"FAILED","HUNANA Analysis")
		fileCleanUp(listOfFiles)
		logging.exception('There was an error at HUNANA')
		return False
		#exit(1)

	jobUpdate(id,"EXECUTING","ATAT Analysis")
	ATAT=Atat()
	success=ATAT.run(id=baseFile,AnimalFastaFile=animalFastaFile,AnimalHunanaData=animalHunanaResult,HumanFastaFile=humanFastaFile,HumanHunanaData=humanHunanaResult)
	if not success:
		jobUpdate(id,"FAILED","ATAT Analysis")
		fileCleanUp(listOfFiles)
		logging.exception('There was en error at ATAT')
		return False
	else:
		fileCleanUp(listOfFiles)
		jobUpdate(id,"DONE","Analysis Complete")
		return True
