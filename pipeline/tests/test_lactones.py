"""Tests for lactone enumeration and reactions."""

import pytest

from pipeline.enumerate.monosaccharides import enumerate_all_monosaccharides
from pipeline.enumerate.sugar_acids import generate_sugar_acids
from pipeline.enumerate.lactones import generate_lactones, CURATED_LACTONES
from pipeline.reactions.lactone_reactions import generate_lactonizations


@pytest.fixture(scope="module")
def monosaccharides():
    return enumerate_all_monosaccharides()


@pytest.fixture(scope="module")
def sugar_acids(monosaccharides):
    return generate_sugar_acids(monosaccharides)


@pytest.fixture(scope="module")
def lactones(sugar_acids):
    return generate_lactones(sugar_acids)


def test_generate_returns_list(lactones):
    assert isinstance(lactones, list)
    assert len(lactones) == len(CURATED_LACTONES)


def test_all_have_lactone_type(lactones):
    for c in lactones:
        assert c["type"] == "lactone"


def test_gdl_formula(lactones):
    """GDL: gluconic acid C6H12O7 - H2O = C6H10O6."""
    gdl = next(c for c in lactones if c["id"] == "D-GDL")
    assert gdl["formula"] == "C6H10O6"


def test_glucuronolactone_formula(lactones):
    """Glucuronolactone: glucuronic acid C6H10O7 - H2O = C6H8O6."""
    glcal = next(c for c in lactones if c["id"] == "D-GlcAL")
    assert glcal["formula"] == "C6H8O6"


def test_unique_ids(lactones):
    ids = [c["id"] for c in lactones]
    assert len(ids) == len(set(ids))


def test_lactonizations_generated(lactones):
    rxns = generate_lactonizations(lactones)
    assert len(rxns) > 0
    # Should have both forward and reverse for each lactone
    assert len(rxns) == len(lactones) * 2


def test_lactonization_reaction_types(lactones):
    rxns = generate_lactonizations(lactones)
    types = {r["reaction_type"] for r in rxns}
    assert "lactonization" in types
    assert "hydrolysis" in types


def test_lactonization_unique_ids(lactones):
    rxns = generate_lactonizations(lactones)
    ids = [r["id"] for r in rxns]
    assert len(ids) == len(set(ids))
