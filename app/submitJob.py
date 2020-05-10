from flask_wtf import FlaskForm
from wtforms import TextField, SelectField, DecimalField, RadioField, SubmitField, FloatField

# Import Form validators
from wtforms.validators import Required

from elasticsearch import Elasticsearch


class SubmitJobForm(FlaskForm):

	def getProteinList():
        	listOfProteins=[]
        	body = { "size": 0, "aggs": { "protein": { "terms": { "field": "geneProductName.keyword" } } } }
        	#elastic_host=["viva-elastic-1","viva-elastic-2","viva-elastic-3"]
        	elastic_host=["192.168.103.57"]
        	es = Elasticsearch(elastic_host,timeout=1000)
        	data = es.search(index='h5n1',body=body)
        	for protein in data['aggregations']['protein']['buckets']:
                	item=(protein['key'],protein['key'])
                	listOfProteins.append(item)
        	return listOfProteins	

	def getHostList():
                listOfHosts=[]
                body = { "size": 0, "aggs": { "host": { "terms": { "field": "host.keyword" } } } }
                #elastic_host=["viva-elastic-1","viva-elastic-2","viva-elastic-3"]
                elastic_host=["192.168.103.57"]
                es = Elasticsearch(elastic_host,timeout=1000)
                data = es.search(index='h5n1',body=body)
                item=('Avian','Avian')
                listOfHosts.append(item)
                for host in data['aggregations']['host']['buckets']:
                        if host['key']=='Human':
                            continue
                        parsedHost=host['key'].split('/')
                        item=(parsedHost[0],parsedHost[0])
                        listOfHosts.append(item)
                return listOfHosts

	index = SelectField('Index',choices=[('h5n1')])
	protein = SelectField('Protein Type',choices=getProteinList())
	year = TextField('Year Data')
	gapThreshold = TextField('Gap Threshold')
	#hostSpecies = TextField('Host Species')
	hostSpecies = SelectField('Host',default='Avian',choices=getHostList())
	removeRedundant = RadioField('Remove Redundant',default='no',choices=[('yes','Yes'),('no','No')])
