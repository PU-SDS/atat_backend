import json
import operator
import logging
from Bio import SeqIO
from operator import itemgetter
from elasticsearch import Elasticsearch
import elasticsearch.helpers as eshelpers

class Atat:

	hostFastaData=[]
	reservoirFastaData=[]
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
		fastaData=self.reservoirFastaData
		if dataSet=='host':
			fastaData=self.hostFastaData
		for data in fastaData:
			if sequences in str(data.seq):
				listOfHeader=str(data.description).split('|')
				for headerItem in listOfHeader:
					item=str(headerItem).split(':')
					if item[0]=='countryCode' and dataSet=='host':
						listOfSource.append(item[1])
						break
					if item[0]=='host' and dataSet=='reservoir':
						listOfSource.append(item[1])
						break
		return listOfSource

	def getStrain(self,sequences,dataSet):
		listOfStrain=[]
		fastaData=self.reservoirFastaData
		if dataSet=='host':
			fastaData=self.hostFastaData
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
                fastaData=self.reservoirFastaData
                if dataSet=='host':
                        fastaData=self.hostFastaData
                for data in fastaData:
                        if sequences in str(data.seq):
                                listOfHeader=str(data.description).split('|')
                                for headerItem in listOfHeader:
                                        item=str(headerItem).split(':')
                                        if item[0]=='id':
                                                listOfId.append(item[1])
                                                break
                return listOfId


	def run(self,id="",reservoirFastaFile="",reservoirHunanaData=[],hostFastaFile="",hostHunanaData=[]):
		index='flua2h_'+str(id)
		reservoirFastaFileStream=open(reservoirFastaFile)
		hostFastaFileStream=open(hostFastaFile)


		self.reservoirFastaData=list(SeqIO.parse(reservoirFastaFileStream, "fasta"))
		reservoirFastaFileStream.close()
		reservoirHunanaProcessedDataList,reservoirTotalData=self.processHunanaData(reservoirHunanaData)
		SimplifiedreservoirHunanaProcessedDataList=self.simplifyHunanaData(reservoirTotalData,reservoirHunanaProcessedDataList)
		reservoirHunanaProcessedDataListWithMotif=self.motifIdentification(SimplifiedreservoirHunanaProcessedDataList)

		self.hostFastaData=list(SeqIO.parse(hostFastaFileStream, "fasta"))
		hostFastaFileStream.close()
		hostHunanaProcessedDataList,hostTotalData=self.processHunanaData(hostHunanaData)
		SimplifiedhostHunanaProcessedDataList=self.simplifyHunanaData(hostTotalData,hostHunanaProcessedDataList)
		hostHunanaProcessedDataListWithMotif=self.motifIdentification(SimplifiedhostHunanaProcessedDataList)

		finalResult=[]
		for reservoirRecord in reservoirHunanaProcessedDataListWithMotif:
			for hostRecord in hostHunanaProcessedDataListWithMotif:
				if hostRecord['position']==reservoirRecord['position']:
					for reservoirSequence in reservoirRecord['sequences']:
						data={"_index": index,"_type": "master",'Position':reservoirRecord['position'],'Sequence':"",'reservoir':{},'host':{}}
						found=False
						for hostSequence in hostRecord['sequences']:
							if reservoirSequence['sequence']==hostSequence['sequence']:
								data['Sequence']=hostSequence['sequence']
								reservoirData = {\
									'Motif':[reservoirSequence['motifShort'],reservoirSequence['motifLong']],\
									'SeqCount':reservoirRecord['totalSupport'],\
									'SupportPercentage':reservoirSequence['supportPercentage'],\
									'Strain':self.getId(reservoirSequence['sequence'],'reservoir'),\
									'Source':self.getSource(reservoirSequence['sequence'],'reservoir')\
									}
								hostData = {\
									'Motif':[hostSequence['motifShort'],hostSequence['motifLong']],\
									'SeqCount':hostRecord['totalSupport'],\
									'SupportPercentage':hostSequence['supportPercentage'],\
									'Strain':self.getId(hostSequence['sequence'],'host'),\
									'Source':self.getSource(hostSequence['sequence'],'host')\
									}
								data['reservoir']=reservoirData
								data['host']=hostData
								found=True
								finalResult.append(data)
								break
						if not found:
							data['Sequence']=reservoirSequence['sequence']
							reservoirData = {\
								'Motif':[reservoirSequence['motifShort'],reservoirSequence['motifLong']],\
								'SeqCount':reservoirRecord['totalSupport'],\
								'SupportPercentage':reservoirSequence['supportPercentage'],\
								'Strain':self.getId(reservoirSequence['sequence'],'reservoir'),\
								'Source':self.getSource(reservoirSequence['sequence'],'reservoir')\
								}
							hostData = {\
								'Motif':['X','None'],\
								'SeqCount':0,\
								'SupportPercentage':0.0,\
								'Strain':[],\
								'Source':[]\
								}
							data['reservoir']=reservoirData
							data['host']=hostData
							finalResult.append(data)

		for hostRecord in hostHunanaProcessedDataListWithMotif:
			for reservoirRecord in reservoirHunanaProcessedDataListWithMotif:
				if hostRecord['position']==reservoirRecord['position']:
					for hostSequence in hostRecord['sequences']:
						data={"_index": index,"_type": "master",'Position':reservoirRecord['position'],'Sequence':"",'reservoir':{},'host':{}}
						found=False
						for reservoirSequence in reservoirRecord['sequences']:
							if reservoirSequence['sequence']==hostSequence['sequence']:
								found=True
								break
						if not found:
							data['Sequence']=hostSequence['sequence']
							hostData = {\
								'Motif':[hostSequence['motifShort'],hostSequence['motifLong']],\
								'SeqCount':hostRecord['totalSupport'],\
								'SupportPercentage':hostSequence['supportPercentage'],\
								'Strain':self.getId(hostSequence['sequence'],'host'),\
								'Source':self.getSource(hostSequence['sequence'],'host')\
								}
							reservoirData = {\
								'Motif':['X','None'],\
								'SeqCount':0,\
								'SupportPercentage':0.0,\
								'Strain':[],\
								'Source':[]\
								}
							data['reservoir']=reservoirData
							data['host']=hostData
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
