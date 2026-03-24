"""Generate sugar lactone derivatives from sugar acids.

Lactones are cyclic esters formed by intramolecular dehydration of sugar acids.
Formula change: acid - H2O = net -2H, -1O.
"""

import re


def _parse_formula(formula: str) -> dict[str, int]:
    atoms: dict[str, int] = {}
    for match in re.finditer(r'([A-Z][a-z]?)(\d*)', formula):
        element = match.group(1)
        count = int(match.group(2)) if match.group(2) else 1
        if element:
            atoms[element] = atoms.get(element, 0) + count
    return atoms


def _format_formula(atoms: dict[str, int]) -> str:
    order = ["C", "H", "N", "O", "P", "S"]
    parts = []
    for elem in order:
        if elem in atoms and atoms[elem] > 0:
            parts.append(f"{elem}{atoms[elem]}" if atoms[elem] > 1 else elem)
    for elem in sorted(atoms):
        if elem not in order and atoms[elem] > 0:
            parts.append(f"{elem}{atoms[elem]}" if atoms[elem] > 1 else elem)
    return "".join(parts)


def _lactone_formula(acid_formula: str) -> str:
    """Lactone = acid - H2O. Net: -2H, -1O."""
    atoms = _parse_formula(acid_formula)
    atoms["H"] = atoms.get("H", 0) - 2
    atoms["O"] = atoms.get("O", 0) - 1
    return _format_formula(atoms)


# Curated lactones: (acid_id, id, name, aliases)
CURATED_LACTONES = [
    ("D-GlcnA", "D-GDL", "D-Glucono-delta-lactone", ["GDL"]),
    ("D-GlcA", "D-GlcAL", "D-Glucuronolactone", []),
    ("D-GalnA", "D-GalL", "D-Galactonolactone", []),
    ("L-GalA", "L-GulL", "L-Gulonolactone", ["vitamin C precursor"]),
]


def generate_lactones(sugar_acids: list[dict]) -> list[dict]:
    """Generate lactone derivatives from sugar acids.

    Lactones are cyclic esters formed by intramolecular dehydration.
    Each lactone is derived from a parent sugar acid.

    Args:
        sugar_acids: list of sugar acid compound dicts

    Returns:
        list of lactone compound dicts
    """
    acid_map = {c["id"]: c for c in sugar_acids}
    lactones: list[dict] = []

    for acid_id, compound_id, name, aliases in CURATED_LACTONES:
        acid = acid_map.get(acid_id)
        if acid is None:
            raise ValueError(f"Lactone parent acid '{acid_id}' not found")

        compound = {
            "id": compound_id,
            "name": name,
            "aliases": aliases,
            "type": "lactone",
            "carbons": acid["carbons"],
            "chirality": acid["chirality"],
            "formula": _lactone_formula(acid["formula"]),
            "stereocenters": list(acid["stereocenters"]),
            "modifications": [{"type": "lactone", "position": 1}],
            "parent_monosaccharide": acid.get("parent_monosaccharide"),
            "commercial": False,
            "cost_usd_per_kg": None,
            "metadata": {
                "parent_acid": acid_id,
                "acid_type": acid.get("metadata", {}).get("acid_type", ""),
            },
            "chebi_id": None,
            "kegg_id": None,
            "pubchem_id": None,
            "inchi": None,
            "smiles": None,
        }
        lactones.append(compound)

    return lactones
