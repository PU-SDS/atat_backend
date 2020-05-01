from flask import Blueprint, jsonify, request, render_template, \
                  flash, g, session, redirect, url_for

# Import module forms
from app.mod_flua2h.submitJob import SubmitJobForm
from app.mod_flua2h.viewResult import getResult,getPosition
from app.mod_flua2h.viewStatus import getStatus
import requests

# Define the blueprint: 'auth', set its url prefix: app.url/auth
mod_flua2h = Blueprint('flua2h', __name__, url_prefix='/flua2h')

# Set the route and accepted methods
@mod_flua2h.route('/submitJob/', methods=['POST'])
def submitJob():
	form = SubmitJobForm(request.form)
	index = form.index.data
	protein = form.protein.data
	year = form.year.data
	hostSpecies = form.hostSpecies.data
	gapThreshold = form.gapThreshold.data
	removeRedundant = form.removeRedundant.data
	payload = {'index':index, 'hostSpecies':hostSpecies, 'protein':protein , 'gapThreshold':gapThreshold, 'removeRedundant':removeRedundant}
	r = requests.post("http://viva-compute:8080/flua2h", data=payload)
	return jsonify(r.text)
	#return render_template('flua2h/parameter.html', protein=protein, year=year, gapthreshold=gapThreshold, removeredundant=removeRedundant)

@mod_flua2h.route('/submitForm', methods=['GET'])
def submitForm():
	form = SubmitJobForm(request.form)
	return render_template('flua2h/submitJob.html',form=form)

@mod_flua2h.route('/showResult', methods=['GET','POST'])
def showResult():
	position=1
	jobId=''
	if request.method == 'POST':
		#input = request.form
		#jobId = input['jobId']
		jobId = request.form.get('jobId',default='',type=str)
		position = request.form.get('position', default=1, type=int)
	else:
		jobId = request.args.get('jobId',default='',type=str)
		position = request.args.get('position', default=1, type=int)
	result = getResult(jobId,position)
	for item in result:
		AnimalSourceDic={}
		HumanSourceDic={}
		animalTotal=0
		humanTotal=0
		for source in item['Animal']['Source']:
			rsource=source.split('/')
			rsource=rsource[0]
			if rsource in AnimalSourceDic:
				AnimalSourceDic[rsource]=AnimalSourceDic[rsource]+1
			else:
				AnimalSourceDic[rsource]=1
			animalTotal=animalTotal+1
		for source in AnimalSourceDic:
			data=[]
			data.append(AnimalSourceDic[source])
			data.append(round(AnimalSourceDic[source]*100/animalTotal,2))
			AnimalSourceDic[source]=data
		for source in item['Human']['Source']:
			if source in HumanSourceDic:
				HumanSourceDic[source]=HumanSourceDic[source]+1
			else:
				HumanSourceDic[source]=1
			humanTotal=humanTotal+1
		for source in HumanSourceDic:
			data=[]
			data.append(HumanSourceDic[source])
			data.append(round(HumanSourceDic[source]*100/humanTotal,2))
			HumanSourceDic[source]=data
		item['Animal']['Source']=AnimalSourceDic
		item['Human']['Source']=HumanSourceDic
	listOfPosition = getPosition(jobId)
	current_position = position
	listOfMotif = ['Index','Major','Minor','Unique']
	return render_template('flua2h/showResult.html',data=result,position=listOfPosition,motif=listOfMotif,current_position=current_position,jobId=jobId)
	#return jsonify(data)

@mod_flua2h.route('/showStatus', methods=['POST'])
def showStatus():
        input = request.form
        data = getStatus(input['jobId'])
        return jsonify(data)

@mod_flua2h.route('/', methods=['GET'])
def index():
	return render_template('flua2h/index.html')


