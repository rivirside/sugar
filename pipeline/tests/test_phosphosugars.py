"""Tests for phosphosugar enumeration and reactions."""

import re
from pipeline.enumerate.monosaccharides import enumerate_all_monosaccharides
from pipeline.enumerate.phosphosugars import generate_phosphosugars


def _get_all_compounds():
    """Helper: enumerate monosaccharides + phosphosugars."""
    mono = enumerate_all_monosaccharides()
    phospho = generate_phosphosugars(mono)
    return mono, phospho


# --- Enumeration tests ---


def test_systematic_count():
    """136 systematic phosphosugars: 64 aldo-mono + 32 aldo-bis + 32 keto-mono + 8 keto-bis."""
    mono, phospho = _get_all_compounds()
    systematic = [p for p in phospho if not p["metadata"].get("curated")]
    assert len(systematic) == 136


def test_curated_count():
    """8 curated phosphosugars (G3P, DHAP, E4P, R5P, Ru5P, Xu5P, S7P, F26BP)."""
    mono, phospho = _get_all_compounds()
    curated = [p for p in phospho if p["metadata"].get("curated")]
    assert len(curated) == 8


def test_total_count():
    """144 total phosphosugars."""
    mono, phospho = _get_all_compounds()
    assert len(phospho) == 144


def test_all_type_phosphate():
    """Every phosphosugar has type='phosphate'."""
    mono, phospho = _get_all_compounds()
    assert all(p["type"] == "phosphate" for p in phospho)


def test_ids_unique():
    """All phosphosugar IDs are unique."""
    mono, phospho = _get_all_compounds()
    ids = [p["id"] for p in phospho]
    assert len(ids) == len(set(ids)), f"Duplicate IDs: {[x for x in ids if ids.count(x) > 1]}"


def test_modifications_non_empty():
    """Every phosphosugar has a non-empty modifications list."""
    mono, phospho = _get_all_compounds()
    for p in phospho:
        assert isinstance(p["modifications"], list)
        assert len(p["modifications"]) > 0
        for mod in p["modifications"]:
            assert mod["type"] == "phosphate"
            assert isinstance(mod["position"], int)
            assert mod["position"] >= 1


def test_parent_exists():
    """Every parent_monosaccharide references an existing compound."""
    mono, phospho = _get_all_compounds()
    mono_ids = {c["id"] for c in mono}
    for p in phospho:
        assert p["parent_monosaccharide"] in mono_ids, (
            f"{p['id']} references missing parent {p['parent_monosaccharide']}"
        )


def test_stereocenters_inherited():
    """Stereocenters match the parent monosaccharide."""
    mono, phospho = _get_all_compounds()
    mono_map = {c["id"]: c for c in mono}
    for p in phospho:
        parent = mono_map[p["parent_monosaccharide"]]
        assert p["stereocenters"] == parent["stereocenters"], (
            f"{p['id']} stereocenters {p['stereocenters']} != parent {parent['stereocenters']}"
        )


def test_chirality_inherited():
    """Chirality matches the parent monosaccharide."""
    mono, phospho = _get_all_compounds()
    mono_map = {c["id"]: c for c in mono}
    for p in phospho:
        parent = mono_map[p["parent_monosaccharide"]]
        assert p["chirality"] == parent["chirality"]


def test_carbons_inherited():
    """Carbon count matches the parent monosaccharide."""
    mono, phospho = _get_all_compounds()
    mono_map = {c["id"]: c for c in mono}
    for p in phospho:
        parent = mono_map[p["parent_monosaccharide"]]
        assert p["carbons"] == parent["carbons"]


def _parse_formula(formula: str) -> dict[str, int]:
    """Parse molecular formula into element counts."""
    atoms: dict[str, int] = {}
    for match in re.finditer(r'([A-Z][a-z]?)(\d*)', formula):
        element = match.group(1)
        count = int(match.group(2)) if match.group(2) else 1
        if element:
            atoms[element] = atoms.get(element, 0) + count
    return atoms


def test_formula_correctness():
    """Phosphosugar formula = parent + n * PO3H (per phosphate group)."""
    mono, phospho = _get_all_compounds()
    mono_map = {c["id"]: c for c in mono}
    for p in phospho:
        parent = mono_map[p["parent_monosaccharide"]]
        parent_atoms = _parse_formula(parent["formula"])
        phospho_atoms = _parse_formula(p["formula"])
        n_phosphates = len(p["modifications"])

        expected = dict(parent_atoms)
        expected["P"] = expected.get("P", 0) + n_phosphates
        expected["O"] = expected.get("O", 0) + 3 * n_phosphates
        expected["H"] = expected.get("H", 0) + 1 * n_phosphates

        assert phospho_atoms == expected, (
            f"{p['id']}: formula {p['formula']} != expected from {parent['formula']} + {n_phosphates}xPO3H"
        )


def test_metadata_has_required_fields():
    """Metadata includes phosphate_positions, parent_type, curated."""
    mono, phospho = _get_all_compounds()
    for p in phospho:
        meta = p["metadata"]
        assert "phosphate_positions" in meta
        assert "parent_type" in meta
        assert meta["parent_type"] in ("aldose", "ketose")
        assert "curated" in meta


def test_c2_excluded_from_systematic():
    """No systematic phosphosugar has a phosphate at position 2."""
    mono, phospho = _get_all_compounds()
    systematic = [p for p in phospho if not p["metadata"].get("curated")]
    for p in systematic:
        positions = [m["position"] for m in p["modifications"]]
        assert 2 not in positions, f"{p['id']} has position 2 phosphate"


def test_known_compounds_present():
    """Key biologically important phosphosugars are generated."""
    mono, phospho = _get_all_compounds()
    ids = {p["id"] for p in phospho}
    expected = {"D-GLC-6P", "D-FRU-6P", "D-FRU-1,6BP", "D-GLYC-3P", "DHA-1P",
                "D-RIB-5P", "D-RBU-5P", "D-XLU-5P", "D-SED-7P", "D-FRU-2,6BP"}
    missing = expected - ids
    assert not missing, f"Missing expected compounds: {missing}"


from pipeline.enumerate.polyols import generate_polyols
from pipeline.reactions.phosphorylation import (
    generate_phosphorylations,
    generate_dephosphorylations,
    generate_mutases,
    generate_phospho_epimerizations,
    generate_phospho_isomerizations,
)


def _get_all_with_reactions():
    """Helper: get compounds + all phospho-reactions."""
    mono = enumerate_all_monosaccharides()
    phospho = generate_phosphosugars(mono)
    all_compounds = mono + phospho
    phos = generate_phosphorylations(phospho)
    dephos = generate_dephosphorylations(phospho)
    mutases = generate_mutases(phospho)
    epi = generate_phospho_epimerizations(phospho)
    iso = generate_phospho_isomerizations(phospho)
    return all_compounds, phospho, phos, dephos, mutases, epi, iso


# --- Phosphorylation tests ---


def test_phosphorylation_count():
    """One phosphorylation per phosphosugar."""
    _, phospho, phos, _, _, _, _ = _get_all_with_reactions()
    assert len(phos) == len(phospho)


def test_phosphorylation_links_parent():
    """Each phosphorylation substrate is the parent monosaccharide."""
    mono = enumerate_all_monosaccharides()
    phospho = generate_phosphosugars(mono)
    mono_map = {c["id"]: c for c in mono}
    phos = generate_phosphorylations(phospho)
    for r in phos:
        sub_id = r["substrates"][0]
        prod_id = r["products"][0]
        prod = next(p for p in phospho if p["id"] == prod_id)
        assert sub_id == prod["parent_monosaccharide"], (
            f"Reaction {r['id']}: substrate {sub_id} != parent {prod['parent_monosaccharide']}"
        )


def test_phosphorylation_fields():
    """Phosphorylation reactions have correct type, cofactor, and tier."""
    _, _, phos, _, _, _, _ = _get_all_with_reactions()
    for r in phos:
        assert r["reaction_type"] == "phosphorylation"
        assert r["cofactor_burden"] == 1.0
        assert r["evidence_tier"] == "hypothetical"


# --- Dephosphorylation tests ---


def test_dephosphorylation_count():
    """One dephosphorylation per phosphosugar."""
    _, phospho, _, dephos, _, _, _ = _get_all_with_reactions()
    assert len(dephos) == len(phospho)


def test_dephosphorylation_reverse_of_phos():
    """Each dephosphorylation is substrate=phosphosugar, product=parent."""
    mono = enumerate_all_monosaccharides()
    phospho = generate_phosphosugars(mono)
    dephos = generate_dephosphorylations(phospho)
    for r in dephos:
        sub_id = r["substrates"][0]
        prod_id = r["products"][0]
        sub = next(p for p in phospho if p["id"] == sub_id)
        assert prod_id == sub["parent_monosaccharide"]


def test_dephosphorylation_fields():
    """Dephosphorylation has cofactor_burden=0.0."""
    _, _, _, dephos, _, _, _ = _get_all_with_reactions()
    for r in dephos:
        assert r["reaction_type"] == "dephosphorylation"
        assert r["cofactor_burden"] == 0.0


# --- Mutase tests ---


def test_mutase_between_same_parent():
    """Mutases connect mono-phosphosugars sharing the same parent."""
    mono = enumerate_all_monosaccharides()
    phospho = generate_phosphosugars(mono)
    mutases = generate_mutases(phospho)
    phospho_map = {p["id"]: p for p in phospho}
    for r in mutases:
        sub = phospho_map[r["substrates"][0]]
        prod = phospho_map[r["products"][0]]
        assert sub["parent_monosaccharide"] == prod["parent_monosaccharide"]
        assert sub["stereocenters"] == prod["stereocenters"]


def test_mutase_different_positions():
    """Mutase substrate and product have different phosphate positions."""
    mono = enumerate_all_monosaccharides()
    phospho = generate_phosphosugars(mono)
    mutases = generate_mutases(phospho)
    phospho_map = {p["id"]: p for p in phospho}
    for r in mutases:
        sub = phospho_map[r["substrates"][0]]
        prod = phospho_map[r["products"][0]]
        assert sub["metadata"]["phosphate_positions"] != prod["metadata"]["phosphate_positions"]


def test_mutase_only_mono_phosphates():
    """Mutases only involve mono-phosphosugars (not bisphosphates)."""
    mono = enumerate_all_monosaccharides()
    phospho = generate_phosphosugars(mono)
    mutases = generate_mutases(phospho)
    phospho_map = {p["id"]: p for p in phospho}
    for r in mutases:
        sub = phospho_map[r["substrates"][0]]
        prod = phospho_map[r["products"][0]]
        assert len(sub["modifications"]) == 1
        assert len(prod["modifications"]) == 1


def test_mutase_count():
    """Expected: 24 parents x C(4,2) pairs x 2 directions = 288 directed."""
    _, _, _, _, mutases, _, _ = _get_all_with_reactions()
    assert len(mutases) == 288


def test_mutase_fields():
    """Mutase fields are correct."""
    _, _, _, _, mutases, _, _ = _get_all_with_reactions()
    for r in mutases:
        assert r["reaction_type"] == "mutase"
        assert r["cofactor_burden"] == 0.0


# --- Phospho-epimerization tests ---


def test_phospho_epi_single_stereocenter_diff():
    """Phospho-epimerizations differ at exactly one stereocenter."""
    mono = enumerate_all_monosaccharides()
    phospho = generate_phosphosugars(mono)
    epi = generate_phospho_epimerizations(phospho)
    phospho_map = {p["id"]: p for p in phospho}
    for r in epi:
        sub = phospho_map[r["substrates"][0]]
        prod = phospho_map[r["products"][0]]
        diffs = sum(1 for a, b in zip(sub["stereocenters"], prod["stereocenters"]) if a != b)
        assert diffs == 1, f"Epimerization {r['id']} differs at {diffs} centers"


def test_phospho_epi_matching_modifications():
    """Phospho-epimerizations have matching phosphate patterns."""
    mono = enumerate_all_monosaccharides()
    phospho = generate_phosphosugars(mono)
    epi = generate_phospho_epimerizations(phospho)
    phospho_map = {p["id"]: p for p in phospho}
    for r in epi:
        sub = phospho_map[r["substrates"][0]]
        prod = phospho_map[r["products"][0]]
        assert sub["metadata"]["phosphate_positions"] == prod["metadata"]["phosphate_positions"]


def test_phospho_epi_type():
    """Phospho-epimerization uses reaction_type='epimerization'."""
    _, _, _, _, _, epi, _ = _get_all_with_reactions()
    for r in epi:
        assert r["reaction_type"] == "epimerization"


# --- Phospho-isomerization tests ---


def test_phospho_iso_aldose_ketose_pairs():
    """Phospho-isomerizations connect aldose-P to ketose-P."""
    mono = enumerate_all_monosaccharides()
    phospho = generate_phosphosugars(mono)
    iso = generate_phospho_isomerizations(phospho)
    phospho_map = {p["id"]: p for p in phospho}
    for r in iso:
        sub = phospho_map[r["substrates"][0]]
        prod = phospho_map[r["products"][0]]
        types = {sub["metadata"]["parent_type"], prod["metadata"]["parent_type"]}
        assert types == {"aldose", "ketose"}, f"Iso {r['id']} types: {types}"


def test_phospho_iso_matching_modifications():
    """Phospho-isomerizations have matching phosphate patterns."""
    mono = enumerate_all_monosaccharides()
    phospho = generate_phosphosugars(mono)
    iso = generate_phospho_isomerizations(phospho)
    phospho_map = {p["id"]: p for p in phospho}
    for r in iso:
        sub = phospho_map[r["substrates"][0]]
        prod = phospho_map[r["products"][0]]
        assert sub["metadata"]["phosphate_positions"] == prod["metadata"]["phosphate_positions"]


def test_phospho_iso_type():
    """Phospho-isomerization uses reaction_type='isomerization'."""
    _, _, _, _, _, _, iso = _get_all_with_reactions()
    for r in iso:
        assert r["reaction_type"] == "isomerization"


def test_phospho_epi_count():
    """Expected: 506 directed phospho-epimerization reactions.

    Systematic aldohexose groups: 6 patterns x 32 pairs x 2 = 384
    Systematic ketohexose groups: 5 patterns x 12 pairs x 2 = 120
    Curated additions: D-RBU-5P <-> D-XLU-5P (2 directed)
    Total: 506
    """
    _, _, _, _, _, epi, _ = _get_all_with_reactions()
    assert len(epi) == 506


def test_phospho_iso_count():
    """Expected: 162 directed phospho-isomerization reactions.

    4 shared mono-P patterns x 16 x 2 = 128
    1 shared bis-P pattern x 16 x 2 = 32
    Curated additions: D-RIB-5P <-> D-RBU-5P (2 directed)
    Total: 162
    """
    _, _, _, _, _, _, iso = _get_all_with_reactions()
    assert len(iso) == 162


# --- Cross-cutting tests ---


def test_all_reaction_ids_unique():
    """All phosphosugar reaction IDs are unique."""
    _, _, phos, dephos, mutases, epi, iso = _get_all_with_reactions()
    all_rxns = phos + dephos + mutases + epi + iso
    ids = [r["id"] for r in all_rxns]
    assert len(ids) == len(set(ids)), f"Duplicate IDs found"


def test_all_reactions_hypothetical():
    """All generated reactions have evidence_tier='hypothetical'."""
    _, _, phos, dephos, mutases, epi, iso = _get_all_with_reactions()
    all_rxns = phos + dephos + mutases + epi + iso
    for r in all_rxns:
        assert r["evidence_tier"] == "hypothetical"


def test_mass_balance_carbons():
    """All reactions conserve carbon count."""
    mono = enumerate_all_monosaccharides()
    phospho = generate_phosphosugars(mono)
    all_compounds = mono + phospho
    compound_map = {c["id"]: c for c in all_compounds}

    _, _, phos, dephos, mutases, epi, iso = _get_all_with_reactions()
    for r in phos + dephos + mutases + epi + iso:
        sub_carbons = sum(compound_map[s]["carbons"] for s in r["substrates"])
        prod_carbons = sum(compound_map[p]["carbons"] for p in r["products"])
        assert sub_carbons == prod_carbons, f"Reaction {r['id']}: {sub_carbons} != {prod_carbons}"
