from elasticsearch import Elasticsearch

#elastic_host=["viva-elastic-1","viva-elastic-2","viva-elastic-3"]
elastic_host=["192.168.103.57"]	
es = Elasticsearch(elastic_host,timeout=1000)


def getStatus(jobId):
	cleanJobId = jobId.lstrip("flua2h_")
	body = {"query":{ "match": { "_id": cleanJobId }}}
	# Init scroll by search
	data = es.search(index='jobqueue',body=body)
	result=data['hits']['hits'][0]['_source']
	return result

