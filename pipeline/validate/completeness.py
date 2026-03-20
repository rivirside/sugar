"""Validate that the compound set is complete (all expected stereoisomers present)."""

# Expected counts: key -> expected number of compounds
# Monosaccharides: (type, carbons)
# Phosphosugars: ("phosphate", parent_type, carbons, phosphate_positions_tuple)
EXPECTED_COUNTS: dict[tuple, int] = {}

# Aldoses C2-C7: number of stereocenters = max(0, carbons - 2)
for _carbons in range(2, 8):
    _n_chiral = max(0, _carbons - 2)
    EXPECTED_COUNTS[("aldose", _carbons)] = 2 ** _n_chiral

# Ketoses C3-C7: number of stereocenters = max(0, carbons - 3)
for _carbons in range(3, 8):
    _n_chiral = max(0, _carbons - 3)
    EXPECTED_COUNTS[("ketose", _carbons)] = 2 ** _n_chiral

# Phosphosugars: C6 aldohexoses (16 stereoisomers per position pattern)
for _pos in [(1,), (3,), (4,), (6,), (1, 6), (3, 6)]:
    EXPECTED_COUNTS[("phosphate", "aldose", 6, _pos)] = 16

# Phosphosugars: C6 ketohexoses (8 stereoisomers per position pattern)
for _pos in [(1,), (3,), (4,), (6,), (1, 6)]:
    EXPECTED_COUNTS[("phosphate", "ketose", 6, _pos)] = 8


def check_completeness(compounds: list[dict]) -> list[str]:
    """Check that all expected stereoisomers are present.

    Returns a list of warning strings describing any missing compounds.
    An empty list means the set is complete.
    """
    counts: dict[tuple, int] = {}
    for c in compounds:
        ctype = c.get("type")
        carbons = c.get("carbons")

        if ctype in ("aldose", "ketose") and carbons is not None:
            key = (ctype, carbons)
            counts[key] = counts.get(key, 0) + 1
        elif ctype == "phosphate" and not c.get("metadata", {}).get("curated"):
            parent_type = c.get("metadata", {}).get("parent_type")
            positions = tuple(sorted(c.get("metadata", {}).get("phosphate_positions", [])))
            if parent_type and carbons is not None and positions:
                key = ("phosphate", parent_type, carbons, positions)
                counts[key] = counts.get(key, 0) + 1

    warnings = []
    for key, expected in sorted(EXPECTED_COUNTS.items(), key=str):
        actual = counts.get(key, 0)
        if actual != expected:
            warnings.append(
                f"Expected {expected} compound(s) for group {key}, found {actual}"
            )

    return warnings
