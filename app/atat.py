import json
import operator
import logging
from Bio import SeqIO
from operator import itemgetter
from elasticsearch import Elasticsearch
import elasticsearch.helpers as eshelpers

class Atat:

	HumanFastaData=[]
	AnimalFastaData=[]
	id=""
	elastic_host=["localhost"]
	es = Elasticsearch(elastic_host,timeout=10000)

	def gendata(self,index,data):
		for record in data:
			yield {
				"_index": index,
				"_type": "master",
				"data": record,
			}

	def process_hits(self,hits=""):
		result=[]
		for item in hits:
			data['position']=item['_source']['position']
			data['variants']=item['_source']['variants']
			data['support']=item['_source']['supports']
			data['entropy']=item['_source']['entropy']
			data['sequences']=item['_source']['sequences']
			result.append(data)
		return result

	def input_download(self,id="",source=""):
		index="hunana_"+source+"_"+id
		body = {"query":{"match_all":{}}}
		# Check index exists
		if not self.es.indices.exists(index=index):
			print("Index " + index + " not exists")
			exit()
		# Init scroll by search
		data = self.es.search(index=index,scroll='2m',body=body)
		# Get the scroll ID
		sid = data['_scroll_id']
		scroll_size = len(data['hits']['hits'])
		# Before scroll, process current batch of hits
		result=self.process_hits(hits=data['hits']['hits'])
		while scroll_size > 0:
			data = es.scroll(scroll_id=sid, scroll='2m')
			# Process current batch of hits
			result=result+self.process_hits(hits=data['hits']['hits'])
			# Update the scroll ID
			sid = data['_scroll_id']
			# Get the number of results that returned in the last scroll
		scroll_size = len(data['hits']['hits'])
		return result

	def processHunanaData(self,hunanaData):
		hunanaList = []
		highestPosition = 0
		position = 0
		for record in hunanaData:
			if record['sequences']=='NOT_ANALYZED':
				continue
			position=record['position']
			#support=record['supports']
			sequence=""
			for seq,seqCount in record['sequences'].items():
				sequence = seq
				support = int(seqCount[0])
				newRecord = {'position': position, 'sequence': sequence, 'support': support}
				hunanaList.append(newRecord)
			if int(position) > int(highestPosition):
				highestPosition=int(position)
		return hunanaList, highestPosition

	def simplifyHunanaData(self,highestPosition,hunanaList):
		simplifiedHunanaList = []
		for position in range(1,highestPosition+1):
			positionDic = {'position': position}
			positionList = []
			totalSupport = 0
			for record in hunanaList:
				if int(record['position'])==position:
					totalSupport+=int(record['support'])
					simplifiedDic={'sequence': record['sequence'], 'support': record['support']}
					positionList.append(simplifiedDic)

			positionDic['totalSupport']=totalSupport
			positionDic['sequences']=positionList
			simplifiedHunanaList.append(positionDic)


		for record in simplifiedHunanaList:
			sortedList = sorted(record['sequences'],key=itemgetter("support"),reverse=True)
			record['sequences']=sortedList

		return simplifiedHunanaList

	def motifIdentification(self,simplifiedHunanaList):
		for record in simplifiedHunanaList:
			count = 1
			for sequence in record['sequences']:
				if int(count)==1:
					if int(sequence['support'])>1:
						sequence['motifShort']="I"
						sequence['motifLong']="Index"
				elif int(count)==2:
					if int(sequence['support'])>1:
						sequence['motifShort']="Ma"
						sequence['motifLong']="Major"
				else:
					sequence['motifShort']="Mi"
					sequence['motifLong']="Minor"
				if int(sequence['support'])==1:
					sequence['motifShort']="U"
					sequence['motifLong']="Unique"
				percent=round(100*float(sequence['support'])/float(record['totalSupport']),2)
				sequence['supportPercentage']=percent
				count+=1

		return simplifiedHunanaList

	def getSource(self,sequences,dataSet):
		listOfSource=[]
		fastaData=self.AnimalFastaData
		if dataSet=='Human':
			fastaData=self.HumanFastaData
		for data in fastaData:
			if sequences in str(data.seq):
				listOfHeader=str(data.description).split('|')
				for headerItem in listOfHeader:
					item=str(headerItem).split(':')
					if item[0]=='countryCode' and dataSet=='Human':
						listOfSource.append(item[1])
						break
					if item[0]=='host' and dataSet=='Animal':
						listOfSource.append(item[1])
						break
		return listOfSource

	def getStrain(self,sequences,dataSet):
		listOfStrain=[]
		fastaData=self.AnimalFastaData
		if dataSet=='Human':
			fastaData=self.HumanFastaData
		for data in fastaData:
			if sequences in str(data.seq):
				listOfHeader=str(data.description).split('|')
				for headerItem in listOfHeader:
					item=str(headerItem).split(':')
					if item[0]=='strain':
						listOfStrain.append(item[1])
						break
		return listOfStrain

	def getId(self,sequences,dataSet):
                listOfId=[]
                fastaData=self.AnimalFastaData
                if dataSet=='Human':
                        fastaData=self.HumanFastaData
                for data in fastaData:
                        if sequences in str(data.seq):
                                listOfHeader=str(data.description).split('|')
                                for headerItem in listOfHeader:
                                        item=str(headerItem).split(':')
                                        if item[0]=='id':
                                                listOfId.append(item[1])
                                                break
                return listOfId


	def run(self,id="",AnimalFastaFile="",AnimalHunanaData=[],HumanFastaFile="",HumanHunanaData=[]):
		index='flua2h_'+str(id)
		AnimalFastaFileStream=open(AnimalFastaFile)
		HumanFastaFileStream=open(HumanFastaFile)


		self.AnimalFastaData=list(SeqIO.parse(AnimalFastaFileStream, "fasta"))
		AnimalFastaFileStream.close()
		AnimalHunanaProcessedDataList,AnimalTotalData=self.processHunanaData(AnimalHunanaData)
		SimplifiedAnimalHunanaProcessedDataList=self.simplifyHunanaData(AnimalTotalData,AnimalHunanaProcessedDataList)
		AnimalHunanaProcessedDataListWithMotif=self.motifIdentification(SimplifiedAnimalHunanaProcessedDataList)

		self.HumanFastaData=list(SeqIO.parse(HumanFastaFileStream, "fasta"))
		HumanFastaFileStream.close()
		HumanHunanaProcessedDataList,HumanTotalData=self.processHunanaData(HumanHunanaData)
		SimplifiedHumanHunanaProcessedDataList=self.simplifyHunanaData(HumanTotalData,HumanHunanaProcessedDataList)
		HumanHunanaProcessedDataListWithMotif=self.motifIdentification(SimplifiedHumanHunanaProcessedDataList)

		finalResult=[]
		for animalRecord in AnimalHunanaProcessedDataListWithMotif:
			for humanRecord in HumanHunanaProcessedDataListWithMotif:
				if humanRecord['position']==animalRecord['position']:
					for animalSequence in animalRecord['sequences']:
						data={"_index": index,"_type": "master",'Position':animalRecord['position'],'Sequence':"",'Animal':{},'Human':{}}
						found=False
						for humanSequence in humanRecord['sequences']:
							if animalSequence['sequence']==humanSequence['sequence']:
								data['Sequence']=humanSequence['sequence']
								animalData = {\
									'Motif':[animalSequence['motifShort'],animalSequence['motifLong']],\
									'SeqCount':animalRecord['totalSupport'],\
									'SupportPercentage':animalSequence['supportPercentage'],\
									'Strain':self.getId(animalSequence['sequence'],'Animal'),\
									'Source':self.getSource(animalSequence['sequence'],'Animal')\
									}
								humanData = {\
									'Motif':[humanSequence['motifShort'],humanSequence['motifLong']],\
									'SeqCount':humanRecord['totalSupport'],\
									'SupportPercentage':humanSequence['supportPercentage'],\
									'Strain':self.getId(humanSequence['sequence'],'Human'),\
									'Source':self.getSource(humanSequence['sequence'],'Human')\
									}
								data['Animal']=animalData
								data['Human']=humanData
								found=True
								finalResult.append(data)
								break
						if not found:
							data['Sequence']=animalSequence['sequence']
							animalData = {\
								'Motif':[animalSequence['motifShort'],animalSequence['motifLong']],\
								'SeqCount':animalRecord['totalSupport'],\
								'SupportPercentage':animalSequence['supportPercentage'],\
								'Strain':self.getId(animalSequence['sequence'],'Animal'),\
								'Source':self.getSource(animalSequence['sequence'],'Animal')\
								}
							humanData = {\
								'Motif':['X','None'],\
								'SeqCount':0,\
								'SupportPercentage':0.0,\
								'Strain':[],\
								'Source':[]\
								}
							data['Animal']=animalData
							data['Human']=humanData
							finalResult.append(data)

		for humanRecord in HumanHunanaProcessedDataListWithMotif:
			for animalRecord in AnimalHunanaProcessedDataListWithMotif:
				if humanRecord['position']==animalRecord['position']:
					for humanSequence in humanRecord['sequences']:
						data={"_index": index,"_type": "master",'Position':animalRecord['position'],'Sequence':"",'Animal':{},'Human':{}}
						found=False
						for animalSequence in animalRecord['sequences']:
							if animalSequence['sequence']==humanSequence['sequence']:
								found=True
								break
						if not found:
							data['Sequence']=humanSequence['sequence']
							humanData = {\
								'Motif':[humanSequence['motifShort'],humanSequence['motifLong']],\
								'SeqCount':humanRecord['totalSupport'],\
								'SupportPercentage':humanSequence['supportPercentage'],\
								'Strain':self.getId(humanSequence['sequence'],'Human'),\
								'Source':self.getSource(humanSequence['sequence'],'Human')\
								}
							animalData = {\
								'Motif':['X','None'],\
								'SeqCount':0,\
								'SupportPercentage':0.0,\
								'Strain':[],\
								'Source':[]\
								}
							data['Animal']=animalData
							data['Human']=humanData
							finalResult.append(data)


		try:
			if self.es.indices.exists(index=index):
				self.es.indices.delete(index=index, ignore=[400, 404])
			for success, info in eshelpers.parallel_bulk(self.es,iter(finalResult)):
				if not success:
					return False
					#return {"Result":"KO","Index":index}
			#return {"Result":"OK","Index":index}
			return True
		except:
			logging.exception('There was an error at when uploading the result')
			self.es.indices.delete(index=index, ignore=[400, 404])
			return False
