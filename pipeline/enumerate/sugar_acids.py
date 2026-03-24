"""Generate sugar acid derivatives from monosaccharides.

Sugar acids are oxidized forms:
- Aldonic acids: C1 aldehyde oxidized to carboxyl (-CHO -> -COOH)
  Formula change: +1O (net, C1 goes from aldehyde to carboxyl)
- Uronic acids: terminal CH2OH oxidized to carboxyl (-CH2OH -> -COOH)
  Formula change: -2H (terminal carbon loses 2H, gains 1O, but also loses 1O...
  actually: CH2OH -> COOH = -2H +1O... no.
  CH2OH (C,2H,OH) -> COOH (C,O,OH) = lose 2H, no change in O.
  Wait: CH2OH has formula CH3O. COOH has formula CHO2.
  Difference: -2H, +1O. So net: -2H +1O)

Actually let's be precise:
- Aldonic acid: aldehyde (CHO) -> carboxyl (COOH).
  CHO -> COOH: adds one O and no change in H. Net: +1O.

- Uronic acid: primary alcohol (CH2OH) -> carboxyl (COOH).
  CH2OH -> COOH: loses 2H, adds 1O. Net: -2H, +1O.
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


def _aldonic_acid_formula(parent_formula: str) -> str:
    """Aldonic acid: CHO -> COOH. Net: +1O."""
    atoms = _parse_formula(parent_formula)
    atoms["O"] = atoms.get("O", 0) + 1
    return _format_formula(atoms)


def _uronic_acid_formula(parent_formula: str) -> str:
    """Uronic acid: CH2OH -> COOH. Net: -2H, +1O."""
    atoms = _parse_formula(parent_formula)
    atoms["H"] = atoms.get("H", 0) - 2
    atoms["O"] = atoms.get("O", 0) + 1
    return _format_formula(atoms)


# Curated sugar acids: (parent_id, acid_type, id, name, aliases)
# acid_type: "aldonic" (C1 oxidized) or "uronic" (C6/terminal oxidized)
CURATED_SUGAR_ACIDS = [
    # Uronic acids (most biologically important)
    ("D-GLC", "uronic", "D-GlcA", "D-Glucuronic acid", ["glucuronate"]),
    ("D-GAL", "uronic", "D-GalA", "D-Galacturonic acid", ["galacturonate"]),
    ("D-MAN", "uronic", "D-ManA", "D-Mannuronic acid", ["mannuronate"]),
    ("L-IDO", "uronic", "L-IdoA", "L-Iduronic acid", ["iduronate"]),
    # Aldonic acids
    ("D-GLC", "aldonic", "D-GlcnA", "D-Gluconic acid", ["gluconate"]),
    ("D-GAL", "aldonic", "D-GalnA", "D-Galactonic acid", ["galactonate"]),
    # L-forms
    ("L-GLC", "uronic", "L-GlcA", "L-Glucuronic acid", []),
    ("L-GAL", "uronic", "L-GalA", "L-Galacturonic acid", []),
]


def generate_sugar_acids(compounds: list[dict]) -> list[dict]:
    """Generate sugar acid derivatives from monosaccharides.

    Args:
        compounds: list of monosaccharide compounds

    Returns:
        list of sugar acid compound dicts
    """
    compound_map = {c["id"]: c for c in compounds}
    acids: list[dict] = []

    for parent_id, acid_type, compound_id, name, aliases in CURATED_SUGAR_ACIDS:
        parent = compound_map.get(parent_id)
        if parent is None:
            raise ValueError(f"Sugar acid parent '{parent_id}' not found")

        if acid_type == "aldonic":
            formula = _aldonic_acid_formula(parent["formula"])
            mod_position = 1
        else:
            formula = _uronic_acid_formula(parent["formula"])
            mod_position = parent["carbons"]

        modifications = [{"type": acid_type, "position": mod_position}]

        compound = {
            "id": compound_id,
            "name": name,
            "aliases": aliases,
            "type": "acid",
            "carbons": parent["carbons"],
            "chirality": parent["chirality"],
            "formula": formula,
            "stereocenters": list(parent["stereocenters"]),
            "modifications": modifications,
            "parent_monosaccharide": parent_id,
            "commercial": False,
            "cost_usd_per_kg": None,
            "metadata": {
                "acid_type": acid_type,
                "oxidized_position": mod_position,
                "parent_type": parent["type"],
            },
            "chebi_id": None,
            "kegg_id": None,
            "pubchem_id": None,
            "inchi": None,
            "smiles": None,
        }
        acids.append(compound)

    return acids
