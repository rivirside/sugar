"""Tests for NDP-sugar enumeration and reactions."""

import pytest

from pipeline.enumerate.monosaccharides import enumerate_all_monosaccharides
from pipeline.enumerate.ndp_sugars import generate_ndp_sugars, CURATED_NDP_SUGARS
from pipeline.reactions.ndp_reactions import generate_ndp_epimerizations


@pytest.fixture(scope="module")
def monosaccharides():
    return enumerate_all_monosaccharides()


@pytest.fixture(scope="module")
def ndp_sugars(monosaccharides):
    return generate_ndp_sugars(monosaccharides)


def test_generate_returns_list(ndp_sugars):
    assert isinstance(ndp_sugars, list)
    assert len(ndp_sugars) == len(CURATED_NDP_SUGARS)


def test_all_have_nucleotide_sugar_type(ndp_sugars):
    for c in ndp_sugars:
        assert c["type"] == "nucleotide_sugar"


def test_udp_glucose_exists(ndp_sugars):
    ids = {c["id"] for c in ndp_sugars}
    assert "UDP-GLC" in ids


def test_gdp_mannose_exists(ndp_sugars):
    ids = {c["id"] for c in ndp_sugars}
    assert "GDP-MAN" in ids


def test_formulas_contain_phosphorus(ndp_sugars):
    for c in ndp_sugars:
        assert "P" in c["formula"], f"{c['id']} formula {c['formula']} missing P"


def test_formulas_contain_nitrogen(ndp_sugars):
    for c in ndp_sugars:
        assert "N" in c["formula"], f"{c['id']} formula {c['formula']} missing N"


def test_unique_ids(ndp_sugars):
    ids = [c["id"] for c in ndp_sugars]
    assert len(ids) == len(set(ids))


def test_ndp_epimerizations_generated(ndp_sugars):
    epis = generate_ndp_epimerizations(ndp_sugars)
    assert isinstance(epis, list)
    assert len(epis) > 0


def test_udp_glc_gal_epimerization(ndp_sugars):
    """UDP-Glucose and UDP-Galactose should epimerize (famous UDP-glucose 4-epimerase)."""
    epis = generate_ndp_epimerizations(ndp_sugars)
    pairs = {(r["substrates"][0], r["products"][0]) for r in epis}
    assert ("UDP-GLC", "UDP-GAL") in pairs or ("UDP-GAL", "UDP-GLC") in pairs


def test_epimerizations_same_nucleotide(ndp_sugars):
    """Epimerizations should only happen within same nucleotide carrier."""
    ndp_map = {c["id"]: c for c in ndp_sugars}
    epis = generate_ndp_epimerizations(ndp_sugars)
    for r in epis:
        sub = ndp_map[r["substrates"][0]]
        prod = ndp_map[r["products"][0]]
        sub_nuc = sub.get("metadata", {}).get("nucleotide")
        prod_nuc = prod.get("metadata", {}).get("nucleotide")
        assert sub_nuc == prod_nuc, f"{r['id']}: {sub_nuc} != {prod_nuc}"


def test_epimerization_unique_ids(ndp_sugars):
    epis = generate_ndp_epimerizations(ndp_sugars)
    ids = [r["id"] for r in epis]
    assert len(ids) == len(set(ids))
