"""Tests for multi-dimensional substrate similarity scoring."""

from pipeline.analyze.similarity import compute_similarity


def _compound(id, type, carbons, stereocenters, modifications=None):
    """Helper: minimal compound dict for similarity tests."""
    return {
        "id": id,
        "type": type,
        "carbons": carbons,
        "stereocenters": stereocenters,
        "modifications": modifications,
    }


# --- Identical compounds ---

def test_identical_compounds():
    """Identical compounds have similarity 1.0 and all distances 0."""
    a = _compound("D-GLC", "aldose", 6, ["R", "S", "R", "R"])
    result = compute_similarity(a, a)
    assert result["overall"] == 1.0
    assert result["stereocenter_distance"] == 0
    assert result["modification_distance"] == 0.0
    assert result["carbon_count_distance"] == 0
    assert result["type_distance"] == 0.0


# --- Stereocenter distance ---

def test_one_stereocenter_different():
    """One stereocenter flip produces stereocenter_distance=1, high overall."""
    a = _compound("D-GLC", "aldose", 6, ["R", "S", "R", "R"])
    b = _compound("D-MAN", "aldose", 6, ["S", "S", "R", "R"])
    result = compute_similarity(a, b)
    assert result["stereocenter_distance"] == 1
    assert result["overall"] > 0.85  # small penalty

def test_all_stereocenters_different():
    """All 4 stereocenters different -> stereocenter_distance=4."""
    a = _compound("D-GLC", "aldose", 6, ["R", "S", "R", "R"])
    b = _compound("L-GLC", "aldose", 6, ["S", "R", "S", "S"])
    result = compute_similarity(a, b)
    assert result["stereocenter_distance"] == 4

def test_zero_stereocenters_no_division_error():
    """Both compounds with empty stereocenters: stereo_norm=0.0, no crash."""
    a = _compound("GLYCO", "aldose", 2, [])
    b = _compound("GLYCO2", "aldose", 2, [])
    result = compute_similarity(a, b)
    assert result["stereocenter_distance"] == 0
    assert result["overall"] == 1.0

def test_different_length_stereocenters():
    """Different stereocenter array lengths (C5 vs C6)."""
    a = _compound("C5", "aldose", 5, ["R", "S", "R"])
    b = _compound("C6", "aldose", 6, ["R", "S", "R", "R"])
    result = compute_similarity(a, b)
    # stereo_distance = at least 1 (length diff), carbon_count_distance = 1
    assert result["stereocenter_distance"] >= 1
    assert result["carbon_count_distance"] == 1


# --- Modification distance ---

def test_same_modifications():
    """Same phosphate at same position: modification_distance=0."""
    mods = [{"type": "phosphate", "position": 6}]
    a = _compound("D-GLC-6P", "phosphate", 6, ["R", "S", "R", "R"], mods)
    b = _compound("D-MAN-6P", "phosphate", 6, ["S", "S", "R", "R"], mods)
    result = compute_similarity(a, b)
    assert result["modification_distance"] == 0.0

def test_different_phosphate_position():
    """Same type, different position: modification_distance=0.4."""
    a = _compound("X-1P", "phosphate", 6, ["R"], [{"type": "phosphate", "position": 1}])
    b = _compound("X-6P", "phosphate", 6, ["R"], [{"type": "phosphate", "position": 6}])
    result = compute_similarity(a, b)
    assert result["modification_distance"] == 0.4

def test_has_vs_no_modifications():
    """One has phosphate, other has none: modification_distance=0.7."""
    a = _compound("D-GLC-6P", "phosphate", 6, ["R"], [{"type": "phosphate", "position": 6}])
    b = _compound("D-GLC", "aldose", 6, ["R"], None)
    result = compute_similarity(a, b)
    assert result["modification_distance"] == 0.7


# --- Carbon count distance ---

def test_same_carbon_count():
    """Same carbons: carbon_count_distance=0."""
    a = _compound("A", "aldose", 6, ["R"])
    b = _compound("B", "aldose", 6, ["S"])
    result = compute_similarity(a, b)
    assert result["carbon_count_distance"] == 0

def test_different_carbon_count():
    """C3 vs C6: carbon_count_distance=3."""
    a = _compound("A", "aldose", 3, ["R"])
    b = _compound("B", "aldose", 6, ["R", "R", "R", "R"])
    result = compute_similarity(a, b)
    assert result["carbon_count_distance"] == 3


# --- Type distance ---

def test_same_type():
    """Same type: type_distance=0."""
    a = _compound("A", "aldose", 6, ["R"])
    b = _compound("B", "aldose", 6, ["S"])
    result = compute_similarity(a, b)
    assert result["type_distance"] == 0.0

def test_aldose_ketose():
    """Aldose vs ketose: type_distance=0.5."""
    a = _compound("A", "aldose", 6, ["R"])
    b = _compound("B", "ketose", 6, ["R"])
    result = compute_similarity(a, b)
    assert result["type_distance"] == 0.5

def test_aldose_polyol():
    """Aldose vs polyol: type_distance=0.7."""
    a = _compound("A", "aldose", 6, ["R"])
    b = _compound("B", "polyol", 6, ["R"])
    result = compute_similarity(a, b)
    assert result["type_distance"] == 0.7

def test_phosphate_type_uses_parent_logic():
    """Phosphate type: similarity should look at parent_type or treat as its parent type.
    Two phosphate compounds with same underlying structure are similar."""
    a = _compound("D-GLC-6P", "phosphate", 6, ["R", "S", "R", "R"],
                  [{"type": "phosphate", "position": 6}])
    b = _compound("D-MAN-6P", "phosphate", 6, ["S", "S", "R", "R"],
                  [{"type": "phosphate", "position": 6}])
    result = compute_similarity(a, b)
    assert result["type_distance"] == 0.0  # both phosphate -> same type


# --- Composite score ---

def test_maximally_different():
    """Maximally different compounds: low overall score."""
    a = _compound("A", "aldose", 6, ["R", "S", "R", "R"])
    b = _compound("B", "polyol", 3, [], [{"type": "phosphate", "position": 1}])
    result = compute_similarity(a, b)
    assert result["overall"] < 0.3  # very different

def test_overall_range():
    """Overall score is always between 0.0 and 1.0."""
    pairs = [
        (_compound("A", "aldose", 6, ["R", "S"]), _compound("B", "polyol", 2, [])),
        (_compound("A", "aldose", 6, ["R"]), _compound("B", "aldose", 6, ["R"])),
    ]
    for a, b in pairs:
        result = compute_similarity(a, b)
        assert 0.0 <= result["overall"] <= 1.0
