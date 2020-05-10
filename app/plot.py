# Import module forms
import matplotlib.pyplot as plt
import seaborn as sns
import requests
import io
import base64
import numpy as np

listOfMotif = ['Index','Major','Minor','Unique','Nonatypes']


def generateLinePlot(data_frame,tag):
	dataset=data_frame[['Total','Entropy','Position']]
	fig, axes = plt.subplots()
	#sns.lineplot(x='Position',y='Entropy',data=dataset,ax=axes)
	#dataset.plot.area(x='Position',y='Entropy',ax=axes)
	axes.set_title('Entropy and Incidence of Total Variants\nfor Each Aligned K-Mers Positions of '+tag+' viral sequences')
	#axes.yaxis.grid(True)
	axes.set_xlabel('')
	axes.set_ylabel('K-mer entropy (bits)')
	axes.set_ylim(0,15)
	axes.set_yticks(np.arange(0, 11, 2))
	plt.stackplot(dataset['Position'],dataset['Entropy'])
	axes2 = axes.twinx()
	sns.lineplot(x='Position',y='Total',data=dataset,ax=axes2,color='r')
	axes2.set_ylim(0,100)
	axes2.set_ylabel('Total variants (%)')
	img = io.BytesIO()
	plt.savefig(img, format='png')
	img.seek(0)
	plot_url = base64.b64encode(img.read()).decode('ascii')
	return plot_url

def generateViolinPlot(data_frame):
	#listOfPlotUrl=[]
	i = 1
	fig = plt.figure()
	fig.subplots_adjust(hspace=0.4,wspace=0.8)
	for motif in listOfMotif:
		dataset=data_frame[motif]
		#fig, axes = plt.subplots(2,3,i)
		axes = fig.add_subplot(2,3,i)
		sns.violinplot(data=dataset)
		#axes.violinplot(dataset=dataset)
		#axes.set_title('Frequency Distribution of the Indicated Proteome Sequence Motifs')
		#axes.yaxis.grid(True)
		axes.set_ylim(0,100)
		axes.set_yticks(np.arange(0, 101, 25))
		axes.set_xlabel('Protein')
		axes.set_ylabel(motif+" (%)")
		i=i+1
	img = io.BytesIO()
	plt.savefig(img, format='svg')
	img.seek(0)
	plot_url = base64.b64encode(img.read()).decode('ascii')
	#listOfPlotUrl.append(plot_url)
	#return listOfPlotUrl
	return plot_url

def generateDynamicPlot(data_frame,tag):
	listOfPlotUrl=[]
	dataset=data_frame[listOfMotif+['Total']]
	#major
	fig, axes = plt.subplots()
	axes.set_yticks(np.arange(0, 101, 25))
	sns.scatterplot(x='Total',y='Index',data=dataset,ax=axes)
	sns.scatterplot(x='Total',y='Total',data=dataset,ax=axes,color='r')
	sns.scatterplot(x='Total',y='Major',data=dataset,ax=axes,color='b')
	img = io.BytesIO()
	plt.savefig(img, format='png')
	img.seek(0)
	plot_url = base64.b64encode(img.read()).decode('ascii')
	listOfPlotUrl.append(plot_url)
	
	#minor and unique
	fig, axes = plt.subplots()
	axes.set_yticks(np.arange(0, 101, 25))
	sns.scatterplot(x='Total',y='Index',data=dataset,ax=axes)
	sns.scatterplot(x='Total',y='Total',data=dataset,ax=axes,color='r')
	sns.scatterplot(x='Total',y='Minor',data=dataset,ax=axes,color='b')
	sns.scatterplot(x='Total',y='Unique',data=dataset,ax=axes,color='g')
	img = io.BytesIO()
	plt.savefig(img, format='png')
	img.seek(0)
	plot_url = base64.b64encode(img.read()).decode('ascii')
	listOfPlotUrl.append(plot_url)

	#nonatypes
	fig, axes = plt.subplots()
	sns.scatterplot(x='Total',y='Nonatypes',data=dataset,ax=axes)
	img = io.BytesIO()
	plt.savefig(img, format='png')
	img.seek(0)
	plot_url = base64.b64encode(img.read()).decode('ascii')
	listOfPlotUrl.append(plot_url)
	return listOfPlotUrl

def generateScatterPlot(data_frame,tag):
	dataset=data_frame[['Total','Entropy']]
	fig, axes = plt.subplots()
	sns.scatterplot(x='Total',y='Entropy',data=dataset,ax=axes)
	axes.set_title('Relationship between entropy and incidence of total variants \n'+tag+' protein k-mer positions')
	#axes.yaxis.grid(True)
	axes.set_xlabel('Total variants (%)')
	axes.set_ylabel('K-mer entropy (bits)')
	axes.set_ylim(0,10)
	img = io.BytesIO()
	plt.savefig(img, format='png')
	img.seek(0)
	plot_url = base64.b64encode(img.read()).decode('ascii')
	return plot_url

