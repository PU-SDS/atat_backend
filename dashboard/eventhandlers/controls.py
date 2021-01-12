import zipfile
import base64
from typing import List, Dict

import pandas as pd

from app.warehousing.mongodb.read import MongoDBRead
from bson.json_util import dumps

from io import BytesIO


def download_results_json_eventhandler(job_id: str) -> dict:
    job_id = job_id.replace('/', '')

    host_results, reservoir_results = MongoDBRead(job_id).get_json()

    host_results_json = dumps(host_results, indent=2)
    reservoir_results_json = dumps(reservoir_results, indent=2)

    byte_stream = _create_zip_file(
        [
            {'filename': 'host.json', 'content': host_results_json},
            {'filename': 'reservoir.json', 'content': reservoir_results_json}
        ]
    )

    encoded_zip_file = base64.b64encode(byte_stream.getvalue()).decode()

    return dict(content=encoded_zip_file, filename=f'ATAT_Job_{job_id}.zip', base64=True)


def download_results_csv_eventhandler(job_id: str) -> dict:
    job_id = job_id.replace('/', '')

    host_results, reservoir_results = MongoDBRead(job_id).get_json()

    host_results_csv = pd.DataFrame(host_results).to_csv()
    reservoir_results_csv = pd.DataFrame(reservoir_results).to_csv()

    byte_stream = _create_zip_file(
        [
            {'filename': 'host.csv', 'content': host_results_csv},
            {'filename': 'reservoir.csv', 'content': reservoir_results_csv}
        ]
    )

    encoded_zip_file = base64.b64encode(byte_stream.getvalue()).decode()

    return dict(content=encoded_zip_file, filename=f'ATAT_Job_{job_id}.zip', base64=True)


def _create_zip_file(file_list: List[Dict]) -> BytesIO:
    byte_stream = BytesIO()
    zip_file = zipfile.ZipFile(byte_stream, 'w', zipfile.ZIP_DEFLATED)

    for file in file_list:
        zip_file.writestr(file.get('filename'), file.get('content'))

    zip_file.close()

    return byte_stream
