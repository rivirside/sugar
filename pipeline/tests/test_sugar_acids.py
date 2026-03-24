"""Tests for sugar acid enumeration and reactions."""

import pytest

from pipeline.enumerate.monosaccharides import enumerate_all_monosaccharides
from pipeline.enumerate.sugar_acids import generate_sugar_acids, CURATED_SUGAR_ACIDS
from pipeline.reactions.acid_reactions import generate_oxidations, generate_acid_epimerizations


@pytest.fixture(scope="module")
def monosaccharides():
    return enumerate_all_monosaccharides()


@pytest.fixture(scope="module")
def sugar_acids(monosaccharides):
    return generate_sugar_acids(monosaccharides)


def test_generate_returns_list(sugar_acids):
    assert isinstance(sugar_acids, list)
    assert len(sugar_acids) == len(CURATED_SUGAR_ACIDS)


def test_all_have_acid_type(sugar_acids):
    for c in sugar_acids:
        assert c["type"] == "acid"


def test_glucuronic_acid_formula(sugar_acids):
    """Uronic acid: C6H12O6 -> C6H10O7 (net -2H +1O)."""
    glca = next(c for c in sugar_acids if c["id"] == "D-GlcA")
    assert glca["formula"] == "C6H10O7"


def test_gluconic_acid_formula(sugar_acids):
    """Aldonic acid: C6H12O6 -> C6H12O7 (net +1O)."""
    glcna = next(c for c in sugar_acids if c["id"] == "D-GlcnA")
    assert glcna["formula"] == "C6H12O7"


def test_unique_ids(sugar_acids):
    ids = [c["id"] for c in sugar_acids]
    assert len(ids) == len(set(ids))


def test_oxidations_generated(sugar_acids):
    oxs = generate_oxidations(sugar_acids)
    assert len(oxs) > 0
    for r in oxs:
        assert r["reaction_type"] == "oxidation"
        assert "NAD+" in r.get("cofactors", [])


def test_acid_epimerizations_generated(sugar_acids):
    epis = generate_acid_epimerizations(sugar_acids)
    assert isinstance(epis, list)
    # Should have at least GlcA <-> ManA or GlcA <-> GalA
    assert len(epis) > 0


def test_acid_epimerization_same_acid_type(sugar_acids):
    epis = generate_acid_epimerizations(sugar_acids)
    acid_map = {c["id"]: c for c in sugar_acids}
    for r in epis:
        sub = acid_map[r["substrates"][0]]
        prod = acid_map[r["products"][0]]
        sub_type = sub.get("metadata", {}).get("acid_type")
        prod_type = prod.get("metadata", {}).get("acid_type")
        assert sub_type == prod_type


def test_oxidation_unique_ids(sugar_acids):
    oxs = generate_oxidations(sugar_acids)
    ids = [r["id"] for r in oxs]
    assert len(ids) == len(set(ids))
