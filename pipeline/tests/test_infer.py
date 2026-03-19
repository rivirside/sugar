"""Tests for D-to-L mirroring inference."""

import pytest
from pipeline.import_.infer import find_mirror_compound, infer_mirrored_reactions


def _make_compounds():
    return [
        {"id": "D-GLC", "name": "D-Glucose", "type": "aldose", "chirality": "D", "stereocenters": ["R", "S", "S", "R"], "carbons": 6},
        {"id": "L-GLC", "name": "L-Glucose", "type": "aldose", "chirality": "L", "stereocenters": ["S", "R", "R", "S"], "carbons": 6},
        {"id": "D-FRU", "name": "D-Fructose", "type": "ketose", "chirality": "D", "stereocenters": ["S", "S", "R"], "carbons": 6},
        {"id": "L-FRU", "name": "L-Fructose", "type": "ketose", "chirality": "L", "stereocenters": ["R", "R", "S"], "carbons": 6},
    ]


def test_find_mirror_d_to_l():
    compounds = _make_compounds()
    mirror = find_mirror_compound("D-GLC", compounds)
    assert mirror is not None
    assert mirror["id"] == "L-GLC"


def test_find_mirror_l_to_d():
    compounds = _make_compounds()
    mirror = find_mirror_compound("L-FRU", compounds)
    assert mirror is not None
    assert mirror["id"] == "D-FRU"


def test_find_mirror_none():
    compounds = [_make_compounds()[0]]
    mirror = find_mirror_compound("D-GLC", compounds)
    assert mirror is None


def test_infer_mirrored_reaction():
    compounds = _make_compounds()
    reactions = [{"id": "RHEA:10001", "substrates": ["D-GLC"], "products": ["D-FRU"], "reaction_type": "isomerization", "evidence_tier": "validated", "evidence_criteria": [{"source": "rhea", "rhea_id": "RHEA:10001"}], "rhea_id": "RHEA:10001", "ec_number": "5.3.1.9", "cofactor_burden": 0.0}]
    existing_reaction_ids = {"RHEA:10001"}
    inferred = infer_mirrored_reactions(reactions, compounds, existing_reaction_ids)
    assert len(inferred) == 1
    assert inferred[0]["substrates"] == ["L-GLC"]
    assert inferred[0]["products"] == ["L-FRU"]
    assert inferred[0]["evidence_tier"] == "inferred"


def test_no_double_inference():
    compounds = _make_compounds()
    reactions = [{"id": "RHEA:10001", "substrates": ["D-GLC"], "products": ["D-FRU"], "reaction_type": "isomerization", "evidence_tier": "validated", "evidence_criteria": [], "rhea_id": "RHEA:10001", "ec_number": "5.3.1.9", "cofactor_burden": 0.0}]
    existing_reaction_ids = {"RHEA:10001", "INFER-RHEA:10001-L"}
    inferred = infer_mirrored_reactions(reactions, compounds, existing_reaction_ids)
    assert len(inferred) == 0
