from elasticsearch import Elasticsearch
from common.sequence import getTemporaryFileName

def jobRegister(dataIndex="",protein="",host="",app=""):
	elastic_host=["viva-elastic-1","viva-elastic-2","viva-elastic-3"]
	es = Elasticsearch(elastic_host,timeout=1000)
	id = baseFile = getTemporaryFileName(ext='')
	jobId=str(id)
	parameters={'index':dataIndex,'id':jobId,'protein':protein,'host':host,'status':'NEW','app':app,'stage':'INITIAL'}
	#try:
	es.index(index='jobqueue',doc_type='master',id=jobId,body=parameters)
	return {'submissionStatus':'OK','jobId':jobId}
	#except:
	#	return {'submissionStatus':'FAILED','jobId':jobId}

def jobUpdate(id="",status="",stage=""):
	elastic_host=["viva-elastic-1","viva-elastic-2","viva-elastic-3"]
	es = Elasticsearch(elastic_host,timeout=1000)
	parameters={'doc': {'status': status, 'stage': stage }}
	es.update('jobqueue',doc_type='master',id=id,body=parameters)

def jobStatus(id=""):
	elastic_host=["viva-elastic-1","viva-elastic-2","viva-elastic-3"]
	es = Elasticsearch(elastic_host,timeout=1000)
	body = {"query":{ "bool": { "must": [ { "match": { "_id": id } } ] } } }
	data = es.search(index='jobqueue',body=body)
	status = {'status':data['hits']['hits'][0]['_source']['status'],'stage':data['hits']['hits'][0]['_source']['stage']}
	return status

