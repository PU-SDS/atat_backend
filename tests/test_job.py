import json
from uuid import uuid4
from os.path import join, dirname

import pytest

from atat_backend.core.models import JobDBModel, Parameters
from atat_backend.core.tasks import dima, atat_standalone, warehousing
from atat_backend.v1.helpers import HelperMethods
from atat_backend.v1.models import CreateVivaJobRequest, Position, CreateJobParameters


@pytest.fixture
def test_host_sequences():
    return open(join(dirname(__file__), 'data/host_sequences.afa'), 'r').read()


@pytest.fixture
def test_reservoir_sequences():
    return open(join(dirname(__file__), 'data/reservoir_sequences.afa'), 'r').read()


@pytest.fixture
def test_host_dima_results():
    json_file = open(join(dirname(__file__), 'data/host_dima_results.json'), 'r')

    return json.load(json_file)


@pytest.fixture
def test_reservoir_dima_results():
    json_file = open(join(dirname(__file__), 'data/reservoir_dima_results.json'), 'r')

    return json.load(json_file)


def test_create_job():
    parameters = Parameters(kmer_length=9, header_format=['Country', 'Host']).save()
    job = JobDBModel(parameters=parameters).save()
    print(f"Created job id: {job.id}.")


def test_standalone_job_steps(test_host_sequences, test_reservoir_sequences):
    job = JobDBModel.objects.first()
    host_dima = dima(test_host_sequences, job.id, job.parameters.kmer_length, job.parameters.header_format)
    reservoir_dima = dima(test_reservoir_sequences, job.id, job.parameters.kmer_length, job.parameters.header_format)
    atat_results = [json.loads(switch.json()) for switch in atat_standalone([host_dima, reservoir_dima], job.id)]
    warehousing([[[host_dima, reservoir_dima], atat_results]], job.id)


def test_viva_job_steps(test_host_dima_results, test_reservoir_dima_results):
    idx = str(uuid4())
    host_positions = [Position(**position) for position in test_host_dima_results.get('results')]
    reservoir_positions = [Position(**position) for position in test_reservoir_dima_results.get('results')]
    parameters = CreateJobParameters(kmer_length=9, header_format=['Country', 'Host'])

    payload = CreateVivaJobRequest(
        id=idx, host_dima_positions=host_positions, reservoir_dima_positions=reservoir_positions, parameters=parameters
    )
    HelperMethods.create_viva_job(payload, True)
