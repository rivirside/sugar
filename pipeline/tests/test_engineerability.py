"""Tests for engineerability score computation."""

from pipeline.analyze.engineerability import compute_score


def test_direct_enzyme_score_zero():
    """Direct enzyme match: engineerability = 0.0."""
    score, components = compute_score(
        coverage_level="direct",
        best_similarity=1.0,
        ec_family_size=50,
        has_pdb=True,
        num_candidates=1,
    )
    assert score == 0.0
    assert components["coverage_level"] == 0.0


def test_no_candidates_score_one():
    """No candidates at all: engineerability = 1.0."""
    score, components = compute_score(
        coverage_level="none",
        best_similarity=0.0,
        ec_family_size=None,
        has_pdb=False,
        num_candidates=0,
    )
    assert score > 0.95  # near-maximum: family_richness default (0.8) prevents exact 1.0
    assert components["coverage_level"] == 1.0
    assert components["best_similarity"] == 1.0
    assert components["family_richness"] == 0.8  # None default
    assert components["structural_data"] == 1.0


def test_layer1_high_similarity_low_score():
    """Layer 1 + high similarity + large family + PDB = low score."""
    score, components = compute_score(
        coverage_level="cross_substrate_l1",
        best_similarity=0.95,
        ec_family_size=60,
        has_pdb=True,
        num_candidates=3,
    )
    assert score < 0.2
    assert components["structural_data"] == 0.0


def test_layer2_medium_score():
    """Layer 2 + medium similarity = medium score."""
    score, components = compute_score(
        coverage_level="cross_substrate_l2",
        best_similarity=0.6,
        ec_family_size=10,
        has_pdb=False,
        num_candidates=1,
    )
    assert 0.3 < score < 0.7


def test_family_only_high_score():
    """Family only (Layer 3): high but not maximum score."""
    score, components = compute_score(
        coverage_level="family_only",
        best_similarity=0.0,
        ec_family_size=5,
        has_pdb=False,
        num_candidates=1,
    )
    assert score > 0.6
    assert score < 1.0


def test_candidate_count_modulates_penalty():
    """More candidates reduce the coverage penalty."""
    score_1, _ = compute_score("cross_substrate_l1", 0.8, 20, True, num_candidates=1)
    score_5, _ = compute_score("cross_substrate_l1", 0.8, 20, True, num_candidates=5)
    assert score_5 < score_1  # more candidates = easier


def test_ec_family_none_defaults():
    """When ec_family_size is None, family_richness defaults to 0.8."""
    _, components = compute_score("cross_substrate_l1", 0.8, None, False, 1)
    assert components["family_richness"] == 0.8


def test_large_family_low_penalty():
    """Family of 50+ enzymes: family_richness = 0.0."""
    _, components = compute_score("cross_substrate_l1", 0.8, 100, True, 1)
    assert components["family_richness"] == 0.0


def test_score_always_in_range():
    """Engineerability score is always between 0.0 and 1.0."""
    test_cases = [
        ("direct", 1.0, 100, True, 10),
        ("none", 0.0, None, False, 0),
        ("cross_substrate_l1", 0.5, 25, True, 3),
        ("cross_substrate_l2", 0.3, 1, False, 1),
        ("family_only", 0.0, 5, False, 1),
    ]
    for level, sim, fam, pdb, count in test_cases:
        score, _ = compute_score(level, sim, fam, pdb, count)
        assert 0.0 <= score <= 1.0, f"Out of range for {level}: {score}"
