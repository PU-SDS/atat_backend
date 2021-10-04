import json
from os.path import join, dirname

import pytest

from viva_atat.core.models import JobDBModel, Parameters
from viva_atat.core.tasks import dima, atat, warehousing


@pytest.fixture
def test_host_sequences():
    return open(join(dirname(__file__), 'data/host_sequences.afa'), 'r').read()


@pytest.fixture
def test_reservoir_sequences():
    return open(join(dirname(__file__), 'data/reservoir_sequences.afa'), 'r').read()


def test_create_job():
    parameters = Parameters(kmer_length=9, header_format=['Country', 'Host']).save()
    job = JobDBModel(parameters=parameters).save()
    print(f"Created job id: {job.id}.")


def test_job_steps(test_host_sequences, test_reservoir_sequences):
    job = JobDBModel.objects.first()
    host_dima = dima(test_host_sequences, job.id, job.parameters.kmer_length, job.parameters.header_format)
    reservoir_dima = dima(test_reservoir_sequences, job.id, job.parameters.kmer_length, job.parameters.header_format)
    atat_results = [json.loads(switch.json()) for switch in atat([host_dima, reservoir_dima], job.id)]
    warehousing([[[host_dima, reservoir_dima], atat_results]], job.id)
