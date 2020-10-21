from os.path import join
from flask import Blueprint, render_template, request
from app.core import Constants, SaveUploads
from app.tasks.run_atat import run_atat

blueprint = Blueprint('root', __name__)


@blueprint.route('/submit_job', methods=['GET', 'POST'])
def submit_job():
    if not request.method == 'POST' or not request.files[Constants.HOST_SEQUENCE] \
            or not request.files[Constants.RESERVOIR_SEQUENCE]:
        return '{status: "Error", reason: "Request method is invalid or ' \
               'host and reservoir files not provided"}'

    # Get the uploaded host and reservoir files
    host_sequence = request.files[Constants.HOST_SEQUENCE]
    reservoir_sequence = request.files[Constants.RESERVOIR_SEQUENCE]

    # Call the save uploads function which returns the job id
    job_id = SaveUploads(host_sequence, reservoir_sequence).save()

    if not job_id:
        return '{status: "Error", reason: "Error while saving upload files"}'

    # Submit the job to the queue
    run_atat.delay(job_id)

    # Return the job id
    return job_id


@blueprint.route('/retrieve_job')
def retrieve_job():
    if not request.method == 'POST' or not request.form['job_id']:
        return '<h1>Error</h1>'


@blueprint.route('/')
def index():
    return render_template('index.html')
