"""Compute engineerability score for reactions.

Combines four components:
1. Coverage level (0.4 weight) - direct, cross-substrate, family, none
2. Best candidate similarity (0.3 weight) - inverted similarity score
3. EC family richness (0.15 weight) - normalized family size
4. Structural data availability (0.15 weight) - PDB crystal structure exists
"""

# Coverage level base penalties
_COVERAGE_PENALTIES = {
    "direct": 0.0,
    "cross_substrate_l1": 0.3,
    "cross_substrate_l2": 0.6,
    "family_only": 0.8,
    "none": 1.0,
}

# Component weights (sum = 1.0)
W_COVERAGE = 0.4
W_SIMILARITY = 0.3
W_FAMILY = 0.15
W_STRUCTURAL = 0.15

# Default family richness penalty when BRENDA data unavailable
DEFAULT_FAMILY_PENALTY = 0.8

# Family size normalization ceiling
FAMILY_SIZE_CEILING = 50


def compute_score(
    coverage_level: str,
    best_similarity: float,
    ec_family_size: int | None,
    has_pdb: bool,
    num_candidates: int = 0,
) -> tuple[float, dict]:
    """Compute the engineerability score for a reaction.

    Args:
        coverage_level: One of "direct", "cross_substrate_l1",
            "cross_substrate_l2", "family_only", "none".
        best_similarity: Overall similarity score of the best candidate
            (0.0 to 1.0). Set to 0.0 when no candidates exist.
        ec_family_size: Number of distinct enzymes in the EC subclass.
            None when BRENDA data is unavailable.
        has_pdb: Whether the best candidate has a PDB crystal structure.
        num_candidates: Number of cross-substrate candidates found.

    Returns:
        Tuple of (engineerability_score, components_dict).
        Score range: 0.0 (trivially engineerable) to 1.0 (no leads).
    """
    # Component 1: Coverage level with candidate count modulation
    base_penalty = _COVERAGE_PENALTIES.get(coverage_level, 1.0)
    if coverage_level != "direct" and coverage_level != "none" and num_candidates > 0:
        modulation = 1.0 - min(num_candidates / 5, 1.0) * 0.3
        coverage_component = base_penalty * modulation
    else:
        coverage_component = base_penalty

    # Component 2: Best candidate similarity (inverted: higher = harder)
    if coverage_level == "direct":
        similarity_component = 0.0
    elif coverage_level == "none":
        similarity_component = 1.0
    else:
        similarity_component = 1.0 - best_similarity

    # Component 3: EC family richness
    if ec_family_size is None:
        family_component = DEFAULT_FAMILY_PENALTY
    else:
        family_component = 1.0 - min(ec_family_size / FAMILY_SIZE_CEILING, 1.0)

    # Component 4: Structural data availability (binary in phase 1)
    structural_component = 0.0 if has_pdb else 1.0

    # Composite score
    score = (
        W_COVERAGE * coverage_component
        + W_SIMILARITY * similarity_component
        + W_FAMILY * family_component
        + W_STRUCTURAL * structural_component
    )
    score = max(0.0, min(1.0, score))

    components = {
        "coverage_level": coverage_component,
        "best_similarity": similarity_component,
        "family_richness": family_component,
        "structural_data": structural_component,
    }

    return score, components
