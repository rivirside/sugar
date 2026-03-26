"""Tests for cross-class bridge reactions."""

import pytest

from pipeline.enumerate.monosaccharides import enumerate_all_monosaccharides
from pipeline.enumerate.polyols import generate_polyols
from pipeline.enumerate.phosphosugars import generate_phosphosugars
from pipeline.enumerate.deoxy_sugars import generate_deoxy_sugars
from pipeline.enumerate.amino_sugars import generate_amino_sugars
from pipeline.enumerate.sugar_acids import generate_sugar_acids
from pipeline.enumerate.lactones import generate_lactones
from pipeline.enumerate.ndp_sugars import generate_ndp_sugars
from pipeline.reactions.bridge_reactions import (
    generate_amination_bridges,
    generate_deoxygenation_bridges,
    generate_ndp_activation_bridges,
)


@pytest.fixture(scope="module")
def all_compounds():
    mono = enumerate_all_monosaccharides()
    polyols = generate_polyols(mono)
    phospho = generate_phosphosugars(mono)
    deoxy = generate_deoxy_sugars(mono)
    amino = generate_amino_sugars(mono)
    acids = generate_sugar_acids(mono)
    lact = generate_lactones(acids)
    ndp = generate_ndp_sugars(mono + deoxy + amino + acids)
    return mono + polyols + phospho + deoxy + amino + acids + lact + ndp


@pytest.fixture(scope="module")
def amino_sugars():
    mono = enumerate_all_monosaccharides()
    return generate_amino_sugars(mono)


@pytest.fixture(scope="module")
def deoxy_sugars():
    mono = enumerate_all_monosaccharides()
    return generate_deoxy_sugars(mono)


@pytest.fixture(scope="module")
def ndp_sugars():
    mono = enumerate_all_monosaccharides()
    deoxy = generate_deoxy_sugars(mono)
    amino = generate_amino_sugars(mono)
    acids = generate_sugar_acids(mono)
    return generate_ndp_sugars(mono + deoxy + amino + acids)


# --- Amination bridges ---

def test_amination_bridges_generated(amino_sugars):
    bridges = generate_amination_bridges(amino_sugars)
    assert len(bridges) > 0


def test_amination_skips_nacetyl(amino_sugars):
    """N-acetyl amino sugars should NOT get direct amination bridges."""
    bridges = generate_amination_bridges(amino_sugars)
    product_ids = {r["products"][0] for r in bridges if r["id"].startswith("AMIN")}
    # No N-acetyl compounds should be direct amination targets
    assert "D-GlcNAc" not in product_ids
    assert "D-GalNAc" not in product_ids


def test_amination_glc_to_glcn(amino_sugars):
    bridges = generate_amination_bridges(amino_sugars)
    pairs = {(r["substrates"][0], r["products"][0]) for r in bridges}
    assert ("D-GLC", "D-GlcN") in pairs


def test_amination_has_glutamine_cofactor(amino_sugars):
    bridges = generate_amination_bridges(amino_sugars)
    fwd = [r for r in bridges if r["id"].startswith("AMIN")]
    for r in fwd:
        assert "glutamine" in r["cofactors"]


def test_amination_has_reverse(amino_sugars):
    bridges = generate_amination_bridges(amino_sugars)
    fwd_pairs = {(r["substrates"][0], r["products"][0])
                 for r in bridges if r["id"].startswith("AMIN")}
    rev_pairs = {(r["substrates"][0], r["products"][0])
                 for r in bridges if r["id"].startswith("DEAMIN")}
    for s, p in fwd_pairs:
        assert (p, s) in rev_pairs


def test_amination_unique_ids(amino_sugars):
    bridges = generate_amination_bridges(amino_sugars)
    ids = [r["id"] for r in bridges]
    assert len(ids) == len(set(ids))


# --- Deoxygenation bridges ---

def test_deoxygenation_bridges_generated(deoxy_sugars):
    bridges = generate_deoxygenation_bridges(deoxy_sugars)
    assert len(bridges) > 0


def test_deoxygenation_gal_to_fuc(deoxy_sugars):
    bridges = generate_deoxygenation_bridges(deoxy_sugars)
    pairs = {(r["substrates"][0], r["products"][0]) for r in bridges}
    assert ("L-GAL", "L-FUC") in pairs


def test_deoxygenation_has_nadph(deoxy_sugars):
    bridges = generate_deoxygenation_bridges(deoxy_sugars)
    fwd = [r for r in bridges if r["id"].startswith("DEOX")]
    for r in fwd:
        assert "NADPH" in r["cofactors"]


def test_deoxygenation_has_reverse(deoxy_sugars):
    bridges = generate_deoxygenation_bridges(deoxy_sugars)
    fwd_count = sum(1 for r in bridges if r["id"].startswith("DEOX"))
    rev_count = sum(1 for r in bridges if r["id"].startswith("REOX"))
    assert fwd_count == rev_count


def test_deoxygenation_unique_ids(deoxy_sugars):
    bridges = generate_deoxygenation_bridges(deoxy_sugars)
    ids = [r["id"] for r in bridges]
    assert len(ids) == len(set(ids))


# --- NDP activation bridges ---

def test_ndp_activation_bridges_generated(ndp_sugars, all_compounds):
    bridges = generate_ndp_activation_bridges(ndp_sugars, all_compounds)
    assert len(bridges) > 0


def test_ndp_glc1p_to_udpglc(ndp_sugars, all_compounds):
    """D-GLC-1P should connect to UDP-GLC."""
    bridges = generate_ndp_activation_bridges(ndp_sugars, all_compounds)
    pairs = {(r["substrates"][0], r["products"][0]) for r in bridges}
    assert ("D-GLC-1P", "UDP-GLC") in pairs


def test_ndp_gal1p_to_udpgal(ndp_sugars, all_compounds):
    bridges = generate_ndp_activation_bridges(ndp_sugars, all_compounds)
    pairs = {(r["substrates"][0], r["products"][0]) for r in bridges}
    assert ("D-GAL-1P", "UDP-GAL") in pairs


def test_ndp_activation_has_ntp_cofactor(ndp_sugars, all_compounds):
    bridges = generate_ndp_activation_bridges(ndp_sugars, all_compounds)
    fwd = [r for r in bridges if r["id"].startswith("NDPACT")]
    for r in fwd:
        assert any(c in r["cofactors"] for c in ["UTP", "GTP", "CTP", "dTTP"])


def test_ndp_activation_has_reverse(ndp_sugars, all_compounds):
    bridges = generate_ndp_activation_bridges(ndp_sugars, all_compounds)
    fwd_count = sum(1 for r in bridges if r["id"].startswith("NDPACT"))
    rev_count = sum(1 for r in bridges if r["id"].startswith("NDPREL"))
    assert fwd_count == rev_count


def test_ndp_activation_unique_ids(ndp_sugars, all_compounds):
    bridges = generate_ndp_activation_bridges(ndp_sugars, all_compounds)
    ids = [r["id"] for r in bridges]
    assert len(ids) == len(set(ids))


# --- Cross-type connectivity ---

def test_all_bridges_connect_different_types(amino_sugars, deoxy_sugars, ndp_sugars, all_compounds):
    """Every bridge reaction should connect two different compound types."""
    compound_map = {c["id"]: c for c in all_compounds}
    all_bridges = (
        generate_amination_bridges(amino_sugars) +
        generate_deoxygenation_bridges(deoxy_sugars) +
        generate_ndp_activation_bridges(ndp_sugars, all_compounds)
    )
    for r in all_bridges:
        sub = compound_map.get(r["substrates"][0])
        prod = compound_map.get(r["products"][0])
        if sub and prod:
            assert sub["type"] != prod["type"], (
                f"{r['id']}: same type {sub['type']}"
            )
