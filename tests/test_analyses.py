from itertools import product
from os.path import join, dirname

import pytest
import json

from viva_atat.core.analysis import Analyses

KMER_LENGTH = 9
HEADER_FORMAT = ["Country", "Host"]


def assert_all_properties(model_results, new_results):
    assert new_results.low_support_count == model_results.get('low_support_count')
    assert new_results.protein_name == model_results.get('protein_name')
    assert new_results.sequence_count == model_results.get('sequence_count')
    assert new_results.support_threshold == model_results.get('support_threshold')

    for new_position, old_position in zip(new_results.results, model_results.get('results')):
        assert new_position.position == old_position.get('position')
        assert new_position.support == old_position.get('support')
        assert new_position.low_support == old_position.get('low_support')
        assert new_position.distinct_variants_count == old_position.get('distinct_variants_count')
        assert round(new_position.distinct_variants_incidence) == round(old_position.get('distinct_variants_incidence'))

        if not new_position.variants or not old_position.get('variants'):
            continue

        for new_variant, old_variant in product(new_position.variants, old_position.get('variants')):
            if new_variant.sequence != old_variant.get('sequence'):
                continue

            assert new_variant.sequence == old_variant.get('sequence')
            assert new_variant.motif_long == old_variant.get('motif_long')
            assert new_variant.motif_short == old_variant.get('motif_short')
            assert round(new_variant.incidence) == round(old_variant.get('incidence'))
            assert new_variant.count == old_variant.get('count')
            assert new_variant.metadata == old_variant.get('metadata')

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


@pytest.fixture
def test_at_results():
    json_file = open(join(dirname(__file__), 'data/at_results.json'), 'r')

    return json.load(json_file)


def test_host_dima(test_host_sequences, test_host_dima_results):
    results = Analyses.dima_analysis(test_host_sequences, KMER_LENGTH, HEADER_FORMAT)
    assert_all_properties(test_host_dima_results, results)


def test_reservoir_dima(test_reservoir_sequences, test_reservoir_dima_results):
    results = Analyses.dima_analysis(test_reservoir_sequences, KMER_LENGTH, HEADER_FORMAT)

    assert_all_properties(test_reservoir_dima_results, results)


def test_at_analysis(test_reservoir_sequences, test_host_sequences, test_at_results):
    host_pos_string = str(Analyses.dima_analysis(test_host_sequences, KMER_LENGTH, HEADER_FORMAT).results)
    reservoir_pos_string = str(Analyses.dima_analysis(test_reservoir_sequences, KMER_LENGTH, HEADER_FORMAT).results)

    host_pos_list = json.loads(host_pos_string)
    reservoir_pos_list = json.loads(reservoir_pos_string)

    ata_results = Analyses.at_analysis(host_pos_list, reservoir_pos_list)

    for model_result_switch, new_result_switch in product(test_at_results, ata_results):
        if new_result_switch.sequence != model_result_switch.get('sequence'):
            continue

        assert new_result_switch.sequence == model_result_switch.get('sequence')
        assert new_result_switch.position == model_result_switch.get('position')
        assert new_result_switch.target.value == model_result_switch.get('target')
        assert new_result_switch.source.value == model_result_switch.get('source')
