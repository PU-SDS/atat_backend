from flask import Flask, render_template, request, Response, jsonify
from flask_celery import make_celery
from flask_bootstrap import Bootstrap
from flask_uploads import UploadSet, configure_uploads, ALL
from app.preprocess import preprocess
from app.hunana import Hunana
from app.atat import Atat
from app.viewResult import getResult
from app.sendNotification import send_email_notification
from app.jobManagement import getAnalysisId,registerAtatAnalysis
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

# ATAT Analysis Task Submission
@celery.task(name='__init__.postJob')
def postJob(analysisId,proteinName, hostInputSequenceFile,reservoirInputSequenceFile,removeDuplicates,gapThreshold,userEmail):
	# Register the Task
	registerAtatAnalysis(app,analysisId,proteinName)

	# Preprocess the input file
	(hostAlignedSequenceFile,reservoirAlignedSequenceFile)=preprocess(app,analysisId=analysisId,hostInputSequenceFile=hostInputSequenceFile,reservoirInputSequenceFile=reservoirInputSequenceFile,removeDuplicates=removeDuplicates,gapThreshold=gapThreshold)

	# Run HUNANA and ATAT analysis
	runAnalysis(analysisId,hostInputSequenceFile=hostAlignedSequenceFile,reservoirInputSequenceFile=reservoirAlignedSequenceFile,entropy=True)

        # Notify user if requested
	if userEmail:
		send_email_notification(app,userEmail,analysisId)

# HUNANA and ATAT Analysis Task submission
def runAnalysis(analysisId,hostInputSequenceFile, reservoirInputSequenceFile,entropy):
    HUNANA=Hunana()

    # Run HUNANA Analysis for Host
    hostHunanaResult=HUNANA.run(id=analysisId,inputfile=inputFile,source="host",entropy=entropy)
    hostHunanaOutputStream=open(app.config['UPLOADS_DEFAULT_DEST']+analysisId+"_host.hunana","w")
    json.dump(hostHunanaResult,hostHunanaOutputStream)
    hostHunanaOutputStream.close()

    # Run HUNANA Analysis for Reservoir
    reservoirHunanaResult=HUNANA.run(id=analysisId,inputfile=inputFile,source="reservoir",entropy=entropy)
    reservoirHunanaOutputStream=open(app.config['UPLOADS_DEFAULT_DEST']+analysisId+"_reservoir.hunana","w")
    json.dump(reservoirHunanaResult,reservoirHunanaOutputStream)
    reservoirHunanaOutputStream.close()

    ATAT=Atat()
    atatResult=ATAT.run(id=analysisId,hostFastaFile=hostInputSequenceFile,hostHunanaData=hostHunanaResult,reservoirFastaFile=reservoirInputSequenceFile,reservoirHunanaData=reservoirHunanaResult)
    atatOutputStream=open(app.config['UPLOADS_DEFAULT_DEST']+analysisId+"_"+source+".atat","w")
    json.dump(atatResult,atatOutputStream)
    atatOutputStream.close()
    return True


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
