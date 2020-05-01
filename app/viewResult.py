from elasticsearch import Elasticsearch
import numpy as np

	
#elastic_host=["viva-elastic-1","viva-elastic-2","viva-elastic-3"]
elastic_host=["192.168.103.57"]
es = Elasticsearch(elastic_host,timeout=1000)


def process_hits(hits=""):
	result=[]
	for item in hits:
		result.append(item['_source'])
	return result

def getResult(index="",position=""):
	body = {"sort" : [ { "Animal.SupportPercentage" : {"order" : "desc"}}, "_score"],"query":{ "match": { "Position": position }}}
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
	resultData=process_hits(hits=data['hits']['hits'])
	while scroll_size > 0:
		data = es.scroll(scroll_id=sid, scroll='2m')
    		# Process current batch of hits
		resultData=resultData+process_hits(hits=data['hits']['hits'])
		# Update the scroll ID
		sid = data['_scroll_id']
		# Get the number of results that returned in the last scroll
		scroll_size = len(data['hits']['hits'])
	return resultData

def getPosition(index="",kmer=9):
	body = {"aggs": { "max_position": {"max": {"field": "Position"}},"min_position": {"min": {"field": "Position"}}},"size":0} 
	if not es.indices.exists(index=index):
		print("Index " + index + " not exists")
		exit()
	data = es.search(index=index,body=body)
	max_position = int(data['aggregations']['max_position']['value'])
	min_position = int(data['aggregations']['min_position']['value'])
	listOfPosition=[]
	while min_position <= max_position:
		position=[]
		positionList=range(min_position,min_position+(kmer-1))
		position.append(np.min(positionList))
		#position.append(int(np.median(positionList)))
		position.append(np.max(positionList))
		listOfPosition.append(position)
		min_position=min_position+1
	return listOfPosition
