import math
import time
import queue
import random
import sys, getopt
#import warnings
from collections import OrderedDict, Counter
from multiprocessing import Queue
from billiard.context import Process
from Bio import SeqIO
from operator import itemgetter
from scipy import stats
from decimal import Decimal
import json
import numpy
import elasticsearch.helpers as eshelpers
from elasticsearch import Elasticsearch
#import matplotlib.pyplot as plt

class HunanaSeq:
	# Initial Variable
	#warnings.filterwarnings('error')
	resultList = []
	entropy_iteration = 100
	enable_entropy = False
	elastic_host=["viva-elastic-1","viva-elastic-2","viva-elastic-3"]
	es = Elasticsearch(elastic_host,timeout=1000)

	def gendata(self,index,data):
		for record in data:
			yield {
				"_index": index,
				"_type": "master",
            			"doc": record,
			}
        
	def processRecord(self,input,maxIteration,kmer,enable_entropy):
		result = []
		position = int(numpy.median(list(range(1,kmer+1))))
		while position <= maxIteration:
			count=0
			positionRecord = {'position':position, 'count':count}
			limit=position+kmer
			unwantedchar = set('-XBJZOU')
			for record in input:
				sequence = record.seq[position:limit]
				if any((c in unwantedchar) for c in sequence):
					continue
				else:
					#if positionRecord.has_key(sequence):
					if str(sequence) in positionRecord:
						positionRecord[str(sequence)]=positionRecord[str(sequence)]+1
					else:
						positionRecord[str(sequence)]=1
					positionRecord['count']=positionRecord['count']+1
			if self.enable_entropy:
				normalized_entropy=self.calcNormalizedEntropy(positionRecord,self.entropy_iteration,positionRecord['count'])
				positionRecord['entropy']=normalized_entropy
			else:
				positionRecord['entropy']=0
			result.append(positionRecord)
			position+=1
		return result

	def run(self,id="",inputfile="",outputfile="",alignmentFormat="fasta",kmer=9,protein="",source=""):
		filestream=open(inputfile)
		alignedSeqData=list(SeqIO.parse(filestream, alignmentFormat))
		filestream.close()
		dataLength = len(str(alignedSeqData[0].seq))
		maxIteration = dataLength - kmer
		iterInit = 0
		result=self.processRecord(alignedSeqData,maxIteration,kmer,"false")
		return self.printOutput(id,result,protein,outputfile,source)


	def calcNormalizedEntropy(self,ListOfData="",iteration="",maxSample=""):
		normalizedEntropy = 0
		totalEntropy = []
		inverse = []
		random.seed()
		for iteration in range(iteration):
			#sample = random.choice(range(maxSample))+1
			sample = random.randint(1,maxSample)
			if sample > maxSample:
				sample=maxSample
			totalSample = 0
			sampleDic = {}
			for sampleID in range(sample):
				sequence = ''
				while True:
					sequence = random.choice(list(ListOfData.keys()))
					if sequence == 'position' or sequence == 'count':
						continue
					else:
						totalSample += 1
						#if sampleDic.has_key(sequence):
						if sequence in sampleDic:
							sampleDic[sequence]=sampleDic[sequence]+1
						else:
							sampleDic[sequence]=1
					break
			tTotalEntropy = 0
			#print sampleDic,"\n\n"
			for record in sampleDic:
				#print "record",record
				entropy = float(sampleDic[record])/float(totalSample) * math.log(float(sampleDic[record])/float(totalSample), 2)*-1
				tTotalEntropy = tTotalEntropy + entropy
			totalEntropy.append(float(round(Decimal(tTotalEntropy),2)))
			inverse.append(float(round(Decimal(1)/Decimal(sample),2)))
		#print("totalEntropy="+str(totalEntropy))
		#print "inverse=",inverse 
		#normalizedEntropy = 0
		try:	
			slope, intercept, r_value, p_value, std_err = stats.linregress(inverse,totalEntropy)
			normalizedEntropy = intercept
			#plt.plot(inverse,totalEntropy)
			#plt.show()
		except:
			pass
		return normalizedEntropy		

	def printOutput(self,id="",resultListSorted="",protein="",outputfile="",source=""):
		index='hunana_'+source+"_"+str(id)
		index_config={"index.mapping.total_fields.limit": 100000}
		#("PROTEIN:NO_OF_VARIANTS:POS:SUPPORT:SEQUENCE:COUNT:PERCENTAGE_CONSERVATION:ENTROPY:SUPPORT_GAP\n")
		result=[]
		for recordDic in resultListSorted:
			recordForTheFollowingPosition={}
			#oRecordDic = OrderedDict(recordDic.items)
			if recordDic:
				position = recordDic['position']+5
				count = recordDic['count']
				entropy = recordDic['entropy']
				if count > 0:
					recordDicSorted = sorted(recordDic.items(), key=itemgetter(1),reverse=True)
					firstRecord=True
					sequences={}
					for key,value in recordDicSorted:
						if key == "position" or key == "count" or key == "entropy":
							continue
						percentage_conservation=float(value)*100/float(count)
						if firstRecord:
							recordForTheFollowingPosition['position']=position
							recordForTheFollowingPosition['variants']=len(recordDic)-4
							recordForTheFollowingPosition['supports']=count
							recordForTheFollowingPosition['entropy']=entropy
							sequences[key]=[str(value),str(percentage_conservation)]	
							#filestream.write(protein+":"+str(len(recordDic)-4)+":"+str(position)+":"+str(count)+":"+key+":"+str(value)+":"+str(percentage_conservation)+"%:"+str(entropy)+":0\n")
							firstRecord=False
						else:
							sequences[key]=[str(value),str(percentage_conservation)]
							#filestream.write(":::"+str(position)+":"+key+":"+str(value)+":"+str(percentage_conservation)+"%:"+":0\n")
					recordForTheFollowingPosition['sequences']=sequences
				else:
					recordForTheFollowingPosition['sequences']="NOT_ANALYZED"
					#filestream.write("::NOT_ANALYSED\n")
				#filestream.write(":==========================\n")
				result.append(recordForTheFollowingPosition)
		'''
		if self.es.indices.exists(index=index):
			self.es.index(index=index,doc_type='master',body=recordForTheFollowingPosition)
		else:
			self.es.index(index=index,doc_type='master',body=recordForTheFollowingPosition)
			self.es.indices.put_settings(index=index,body=index_config)
		'''
		#self.es.index(index=index,doc_type='master',body=result)
		#self.es.indices.create(index=index,body=index_config)
		#for success, info in eshelpers.parallel_bulk(self.es,self.gendata(index,result)):
		#	if not success:
		#		return {"Result":"KO","Index":index}
		#return {"Result":"OK","Index":index}
		return result
		#json.dump(result, filestream)
		#filestream.close()
