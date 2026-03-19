"""Tests for ChEBI importer."""

import pytest
from pipeline.import_.chebi import parse_chebi_compounds_tsv, parse_chebi_names_tsv, build_chebi_index


SAMPLE_COMPOUNDS_TSV = """ID\tSTATUS\tSOURCE\tPARENT_ID\tNAME\tDEFINITION\tMODIFIED_ON\tCREATED_BY\tSTAR
17634\tC\tKEGG COMPOUND\t\tD-glucopyranose\tA glucopyranose having D-configuration.\t2023-10-01\tCHEBI\t3
28729\tC\tKEGG COMPOUND\t\tD-mannopyranose\tA mannopyranose having D-configuration.\t2023-10-01\tCHEBI\t3"""

SAMPLE_NAMES_TSV = """ID\tCOMPOUND_ID\tNAME\tTYPE\tSOURCE\tADAPTED\tLANGUAGE
1\t17634\tD-Glucose\tSYNONYM\tChEBI\tfalse\ten
2\t17634\tDextrose\tSYNONYM\tChEBI\tfalse\ten
3\t17634\tGrape sugar\tSYNONYM\tChEBI\tfalse\ten
4\t28729\tD-Mannose\tSYNONYM\tChEBI\tfalse\ten"""


def test_parse_compounds_tsv():
    entries = parse_chebi_compounds_tsv(SAMPLE_COMPOUNDS_TSV)
    assert len(entries) == 2
    assert entries["17634"]["name"] == "D-glucopyranose"
    assert entries["28729"]["name"] == "D-mannopyranose"


def test_parse_names_tsv():
    names = parse_chebi_names_tsv(SAMPLE_NAMES_TSV)
    assert "17634" in names
    assert "D-Glucose" in names["17634"]
    assert "Dextrose" in names["17634"]
    assert "28729" in names
    assert "D-Mannose" in names["28729"]


def test_build_index():
    compounds = parse_chebi_compounds_tsv(SAMPLE_COMPOUNDS_TSV)
    names = parse_chebi_names_tsv(SAMPLE_NAMES_TSV)
    index = build_chebi_index(compounds, names)
    assert "d-glucose" in index
    assert "dextrose" in index
    assert "d-mannose" in index
    assert "d-glucopyranose" in index
    assert index["d-glucose"]["chebi_id"] == "CHEBI:17634"
    assert index["dextrose"]["chebi_id"] == "CHEBI:17634"
