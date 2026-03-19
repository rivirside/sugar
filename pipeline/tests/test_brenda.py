"""Tests for BRENDA importer."""

import pytest
from pipeline.import_.brenda import parse_brenda_km_data, parse_brenda_kcat_data


SAMPLE_KM_RESPONSE = [
    {"ecNumber": "5.3.1.9", "kmValue": 0.5, "substrate": "D-glucose", "organism": "Homo sapiens"},
    {"ecNumber": "5.3.1.9", "kmValue": 1.2, "substrate": "D-glucose", "organism": "Escherichia coli"},
]

SAMPLE_KCAT_RESPONSE = [
    {"ecNumber": "5.3.1.9", "turnoverNumber": 500, "substrate": "D-glucose", "organism": "Homo sapiens"},
]


def test_parse_km_data():
    result = parse_brenda_km_data(SAMPLE_KM_RESPONSE)
    assert len(result) == 2
    assert result[0]["km_mm"] == 0.5
    assert result[0]["organism"] == "Homo sapiens"
    assert result[0]["ec_number"] == "5.3.1.9"


def test_parse_kcat_data():
    result = parse_brenda_kcat_data(SAMPLE_KCAT_RESPONSE)
    assert len(result) == 1
    assert result[0]["kcat_sec"] == 500
    assert result[0]["organism"] == "Homo sapiens"
