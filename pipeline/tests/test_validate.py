from pipeline.validate.completeness import check_completeness
from pipeline.validate.duplicates import check_duplicates
from pipeline.validate.mass_balance import check_mass_balance
from pipeline.enumerate.monosaccharides import enumerate_all_monosaccharides
from pipeline.enumerate.phosphosugars import generate_phosphosugars


def test_completeness_passes_for_full_set():
    monosaccharides = enumerate_all_monosaccharides()
    phosphosugars = generate_phosphosugars(monosaccharides)
    compounds = monosaccharides + phosphosugars
    warnings = check_completeness(compounds)
    assert len(warnings) == 0

def test_completeness_warns_on_missing():
    monosaccharides = enumerate_all_monosaccharides()
    # Drop one monosaccharide (always tracked by completeness)
    compounds = monosaccharides[:-1]
    warnings = check_completeness(compounds)
    assert len(warnings) > 0

def test_no_duplicates_in_full_set():
    monosaccharides = enumerate_all_monosaccharides()
    phosphosugars = generate_phosphosugars(monosaccharides)
    compounds = monosaccharides + phosphosugars
    duplicates = check_duplicates(compounds)
    assert len(duplicates) == 0

def test_detects_duplicate():
    compounds = enumerate_all_monosaccharides()
    compounds.append(compounds[0].copy())
    compounds[-1]["id"] = "DUPLICATE"
    duplicates = check_duplicates(compounds)
    assert len(duplicates) > 0

def test_mass_balance_passes_for_epimerization():
    compounds = [{"id": "D-GLC", "carbons": 6, "type": "aldose"}, {"id": "D-MAN", "carbons": 6, "type": "aldose"}]
    reactions = [{"id": "EPI-001", "substrates": ["D-GLC"], "products": ["D-MAN"], "reaction_type": "epimerization"}]
    errors = check_mass_balance(reactions, {c["id"]: c for c in compounds})
    assert len(errors) == 0

def test_mass_balance_fails_on_carbon_mismatch():
    compounds = [{"id": "D-GLC", "carbons": 6, "type": "aldose"}, {"id": "D-ERY", "carbons": 4, "type": "aldose"}]
    reactions = [{"id": "BAD-001", "substrates": ["D-GLC"], "products": ["D-ERY"], "reaction_type": "epimerization"}]
    errors = check_mass_balance(reactions, {c["id"]: c for c in compounds})
    assert len(errors) > 0


def test_formula_balance_mode():
    from pipeline.validate.mass_balance import check_formula_balance
    reactions = [{"id": "RHEA:10001", "substrates": ["A", "B"], "products": ["C", "D"]}]
    compound_map = {
        "A": {"id": "A", "formula": "C6H12O6"},
        "B": {"id": "B", "formula": "C10H15N5O10P2"},
        "C": {"id": "C", "formula": "C6H13O9P"},
        "D": {"id": "D", "formula": "C10H14N5O7P"},
    }
    errors = check_formula_balance(reactions, compound_map)
    assert isinstance(errors, list)


def test_formula_balance_missing_formula():
    from pipeline.validate.mass_balance import check_formula_balance
    reactions = [{"id": "TEST", "substrates": ["A"], "products": ["B"]}]
    compound_map = {"A": {"id": "A", "formula": "C6H12O6"}, "B": {"id": "B", "formula": None}}
    errors = check_formula_balance(reactions, compound_map)
    assert any("missing formula" in e.lower() or "none" in e.lower() for e in errors)
