from elasticsearch import Elasticsearch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from io import StringIO
from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
import random
import string
import uuid
import os

elastic_host=["viva-elastic-1","viva-elastic-2","viva-elastic-3"]

def getTemporaryFileName(ext=".fasta"):
	random = str(uuid.uuid4())
	random = random.replace("-","")
	return random[0:10]+ext

def process_hits(hits="",speciesType=""):
	sequences=[]
	for item in hits:
		id=item['_id']
		strainName=item['_source']['strainName']
		sequence=item['_source']['aaSeq']
		hostSpecies=item['_source']['host']
		countryCode=item['_source']['gbCountryCode']
		header="strain:"+strainName+"|"+"host:"+hostSpecies+"|"+"countryCode:"+countryCode+"|"+"speciesType:"+speciesType+"|"+"id:"+id
		#record = {'strainName':header,'hostSpecies':hostspecies,'sequence':sequence}
		record = SeqRecord(Seq(sequence,IUPAC.protein),id=header,description="")
		sequences.append(record)
	return sequences

def sequence_download(index="",proteinSymbol="",host="",speciesType=""):
	body = {"query":{ "bool": { "must": [ { "match": { "host": host }}, { "match": { "geneProductName": proteinSymbol }} ] } } }
	es = Elasticsearch(elastic_host,timeout=1000)
	# Check index exists
	if not es.indices.exists(index=index):
		print("Index " + index + " not exists")
		exit()
	# Init scroll by search
	data = es.search(index=index,scroll='2m',body=body)
	# Get the scroll ID
	sid = data['_scroll_id']
	scroll_size = len(data['hits']['hits'])
	# Before scroll, process current batch of hits
	sequences=process_hits(hits=data['hits']['hits'],speciesType=speciesType)
	while scroll_size > 0:
		data = es.scroll(scroll_id=sid, scroll='2m')
    		# Process current batch of hits
		sequences=sequences+process_hits(hits=data['hits']['hits'],speciesType=speciesType)
		# Update the scroll ID
		sid = data['_scroll_id']
		# Get the number of results that returned in the last scroll
		scroll_size = len(data['hits']['hits'])
	return sequences

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
	mafft_cline = MafftCommandline('bioapps/mafft-linux64/mafft.bat',input=file)
	stdout, stderr = mafft_cline()
	align = AlignIO.read(StringIO(stdout), "fasta")
	return align

def sequence_alignmentqc(infile='',outfile='',gap='0.05'):
	status = os.system('bioapps/trimal -keepheader -in '+infile+' -out '+outfile+' -gt '+gap)
	return status

def sequence_split(infile='',prefix='',host='',baseFolder=''):
	os.system("csplit -z "+infile+" -f "+baseFolder+prefix+" '/Human/'")
	os.rename(baseFolder+prefix+"00",baseFolder+prefix+"_Animal.fasta")
	os.rename(baseFolder+prefix+"01",baseFolder+prefix+"_Human.fasta")  
