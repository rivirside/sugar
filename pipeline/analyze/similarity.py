"""Multi-dimensional substrate similarity scoring.

Compares two compound dicts across four dimensions:
1. Stereocenter distance (number of differing positions)
2. Modification distance (phosphate pattern differences)
3. Carbon count distance (chain length difference)
4. Type distance (aldose/ketose/polyol/phosphate)

Returns a dict with per-dimension raw values and a weighted composite score.
"""

# Weights for composite score (sum = 1.0)
W_STEREO = 0.35
W_MOD = 0.25
W_CARBON = 0.25
W_TYPE = 0.15

# Type distance lookup: (type_a, type_b) -> distance
# Symmetric: order doesn't matter.
_TYPE_DISTANCES = {
    ("aldose", "ketose"): 0.5,
    ("aldose", "polyol"): 0.7,
    ("ketose", "polyol"): 0.7,
}


def _stereo_distance(stereo_a: list[str], stereo_b: list[str]) -> int:
    """Count differing stereocenters. Length difference adds to distance."""
    min_len = min(len(stereo_a), len(stereo_b))
    diff = sum(1 for i in range(min_len) if stereo_a[i] != stereo_b[i])
    diff += abs(len(stereo_a) - len(stereo_b))
    return diff


def _modification_distance(mods_a: list[dict] | None, mods_b: list[dict] | None) -> float:
    """Score modification differences.

    Returns:
        0.0: same modifications (identical positions and types)
        0.4: same type, different position
        0.7: different number of modifications (has vs doesn't)
        0.9: different modification types
    """
    a_set = set()
    b_set = set()
    a_types = set()
    b_types = set()

    if mods_a:
        for m in mods_a:
            a_set.add((m["type"], m["position"]))
            a_types.add(m["type"])
    if mods_b:
        for m in mods_b:
            b_set.add((m["type"], m["position"]))
            b_types.add(m["type"])

    # No modifications on either side
    if not a_set and not b_set:
        return 0.0

    # Identical modifications
    if a_set == b_set:
        return 0.0

    # One has mods, other doesn't
    if (not a_set) != (not b_set):
        return 0.7

    # Both have mods but different types
    if a_types != b_types:
        return 0.9

    # Same types, different positions
    return 0.4


def _type_distance(type_a: str, type_b: str) -> float:
    """Score compound type differences."""
    # Normalize: treat "phosphate" as its own type for comparison.
    # Two phosphate compounds = same type (0.0).
    if type_a == type_b:
        return 0.0

    pair = tuple(sorted([type_a, type_b]))
    return _TYPE_DISTANCES.get(pair, 0.5)  # default 0.5 for unknown pairs


def compute_similarity(compound_a: dict, compound_b: dict) -> dict:
    """Compute multi-dimensional substrate similarity between two compounds.

    Args:
        compound_a: First compound dict (must have id, type, carbons,
                     stereocenters, modifications).
        compound_b: Second compound dict.

    Returns:
        Dict with keys:
            overall: float (0.0 = maximally different, 1.0 = identical)
            stereocenter_distance: int (raw count of differing positions)
            modification_distance: float (0.0-1.0)
            carbon_count_distance: int (absolute difference)
            type_distance: float (0.0-1.0)
    """
    stereo_a = compound_a.get("stereocenters", [])
    stereo_b = compound_b.get("stereocenters", [])
    stereo_dist = _stereo_distance(stereo_a, stereo_b)

    max_stereo = max(len(stereo_a), len(stereo_b))
    stereo_norm = stereo_dist / max_stereo if max_stereo > 0 else 0.0

    mod_dist = _modification_distance(
        compound_a.get("modifications"),
        compound_b.get("modifications"),
    )

    carbons_a = compound_a.get("carbons", 0)
    carbons_b = compound_b.get("carbons", 0)
    carbon_dist = abs(carbons_a - carbons_b)
    max_carbons = max(carbons_a, carbons_b)
    carbon_norm = carbon_dist / max_carbons if max_carbons > 0 else 0.0

    type_dist = _type_distance(
        compound_a.get("type", ""),
        compound_b.get("type", ""),
    )

    overall = 1.0 - (
        W_STEREO * stereo_norm
        + W_MOD * mod_dist
        + W_CARBON * carbon_norm
        + W_TYPE * type_dist
    )
    # Clamp to [0.0, 1.0]
    overall = max(0.0, min(1.0, overall))

    return {
        "overall": overall,
        "stereocenter_distance": stereo_dist,
        "modification_distance": mod_dist,
        "carbon_count_distance": carbon_dist,
        "type_distance": type_dist,
    }
