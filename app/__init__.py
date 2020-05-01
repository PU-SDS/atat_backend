from flask import Flask, render_template, request, Response, jsonify
from flask_celery import make_celery
from flask_bootstrap import Bootstrap
from flask_uploads import UploadSet, configure_uploads, ALL
from app.preprocess import preprocess
from app.hunana import Hunana
from app.atat import Atat
from app.viewResult import getResult
from app.sendNotification import send_email_notification
from app.jobManagement import getAnalysisId,registerProteinAnalysis, registerProteomeAnalysis
from app.plot import generateDynamicPlot, generateLinePlot, generateViolinPlot, generateScatterPlot
import json
import time

# Application Initialization
app = Flask(__name__)
app.config.from_object('config')
Bootstrap(app)
celery = make_celery(app)
hostInputSequence = UploadSet('hostInputSequence',ALL)
reservoirInputSequence = UploadSet('reservoirInputSequence',ALL)

configure_uploads(app, (hostInputSequence,reservoirInputSequence))

# Analysis Submission
@app.route('/submit', methods=['GET','POST'])
def submit():
	if request.method == 'POST' and request.files['hostInputSequence'] and request.files['reservoirInputSequence']:
        		# Get Task identifier
		analysisId=getAnalysisId()

		# Upload the input files
		hostInputSequenceFile=hostInputSequence.save(request.files['hostInputSequence'])
        reservoirInputSequenceFile=reservoirInputSequence.save(request.files['reservoirInputSequence'])

		# Set Task parameters
		proteinName=request.form["proteinName"]
		removeDuplicates=request.form["removeDuplicates"]
		gapThreshold=request.form["gapThreshold"]
		if gapThreshold == "":
			gapThreshold = 0.95
		userEmail=request.form["userEmail"]

		# Submit Task to the scheduler
		status=postJob.delay(analysisId,proteinName,hostInputSequenceFile,reservoirInputSequenceFile,removeDuplicates,gapThreshold,userEmail)
		if status:
			message="Your analysis request is successfully submitted , and the analysis id is "+analysisId+". If you provided your email address, you will be notified via email when the analysis is completed."
			return jsonify({'message' : message})
		else:
			message="It seems there was a problem with your analysis request submission. Please check again your input file and/or analysis parameters. Should you need assistance, please do not hesitate to contact us."
			return jsonify({'error' : message})
	else:
		message="It seems there was a problem with your analysis request submission. Please check again your input file and/or analysis parameters. Should you need assistance, please do not hesitate to contact us."
		return jsonify({'error' : message})

# Protein Analysis Task Submission
@celery.task(name='__init__.postJob')
def postJob(analysisId,proteinName, hostInputSequenceFile,reservoirInputSequenceFile,removeDuplicates,gapThreshold,userEmail):
	# Register the Task
	registerAtatAnalysis(app,analysisId,proteinName)

	# Preprocess the input file
	preprocess(app,analysisId=analysisId,inputSequenceFile=inputSequenceFile,removeDuplicates=removeDuplicates,gapThreshold=gapThreshold)

	# HUNANA Analysis
	HUNANA=Hunana()
	inputHunanaResult=HUNANA.run(id=analysisId,inputfile=app.config['UPLOADS_DEFAULT_DEST']+analysisId+"_Input.fasta",source="input",entropy=True)
	outputstream=open(app.config['UPLOADS_DEFAULT_DEST']+analysisId+"_input.hunana","w")
	json.dump(inputHunanaResult,outputstream)
	outputstream.close()

	# VADA Analysis
	ATAT=Atat()
	success=ATAT.run(id=analysisId,InputHunanaData=inputHunanaResult)
	time.sleep(2)

        # Notify user if requested
	if userEmail:
		send_email_notification(app,userEmail,analysisId)

# Result Visualization
@app.route('/viewresult', methods=['GET','POST'])
def viewresult():
	if request.method == 'GET':
		analysisId=request.args.get("analysisId")
	else:
		analysisId=request.form["analysisId"]
	resultInput = getResult("vada_"+analysisId,'Input')
	hosts=['Input']
	return render_template('frontend/result.html',input_violin_plot_url=generateViolinPlot(resultInput),analysisId=analysisId,hosts=hosts)

# Get Raw Result
@app.route('/result/<analysisId>/<host>', methods=['GET','POST'])
def getResultData(analysisId,host):
        result = getResult("vada_"+analysisId,host)
        result = result[['Total','Entropy','Position','Index','Major','Unique','Nonatypes']]
        #return jsonify(result.to_dict())
        #return Response(result.to_csv())
        return Response(result.to_json(orient='records'),mimetype="text/json")


# Main Application Interface
@app.route('/')
def index():
	return render_template('frontend/index.html')

if __name__ == '__main__':
	app.run(debug=True)
