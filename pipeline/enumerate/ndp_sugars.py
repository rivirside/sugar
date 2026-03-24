"""Generate nucleotide sugar (NDP-sugar) derivatives.

NDP-sugars are activated sugar forms where a sugar-1-phosphate is
conjugated to a nucleoside diphosphate (UDP, GDP, CMP, dTDP).

Formula: sugar-1-P + NTP -> NDP-sugar + PPi
But we model NDP-sugars as standalone compounds with their full formula.

Common nucleotide moieties and their formulas (as added to sugar-1-P):
- UDP (uridine diphosphate): adds C9H11N2O8P (UMP minus H2O for condensation)
- GDP (guanosine diphosphate): adds C10H12N5O7P
- CMP (cytidine monophosphate): adds C9H12N3O5 (no extra phosphate)
- dTDP (thymidine diphosphate): adds C10H13N2O8P

Simplified: we store the full molecular formula.
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


# Known NDP-sugar formulas (from biochemistry databases)
# These are curated with exact molecular formulas.
CURATED_NDP_SUGARS = [
    # (id, name, aliases, nucleotide, parent_sugar_id, formula, stereocenters)
    ("UDP-GLC", "UDP-Glucose", ["UDP-D-glucose", "uridine diphosphate glucose"],
     "UDP", "D-GLC", "C15H24N2O17P2", ["R", "S", "R", "R"]),
    ("UDP-GAL", "UDP-Galactose", ["UDP-D-galactose"],
     "UDP", "D-GAL", "C15H24N2O17P2", ["R", "S", "S", "R"]),
    ("UDP-GlcA-NDP", "UDP-Glucuronic acid", ["UDP-glucuronate"],
     "UDP", "D-GlcA", "C15H22N2O18P2", ["R", "S", "R", "R"]),
    ("UDP-GlcNAc-NDP", "UDP-N-Acetylglucosamine", ["UDP-GlcNAc"],
     "UDP", "D-GlcNAc", "C17H27N3O17P2", ["R", "S", "R", "R"]),
    ("UDP-GalNAc-NDP", "UDP-N-Acetylgalactosamine", ["UDP-GalNAc"],
     "UDP", "D-GalNAc", "C17H27N3O17P2", ["R", "S", "S", "R"]),
    ("GDP-MAN", "GDP-Mannose", ["GDP-D-mannose"],
     "GDP", "D-MAN", "C16H25N5O16P2", ["S", "S", "R", "R"]),
    ("GDP-FUC", "GDP-Fucose", ["GDP-L-fucose"],
     "GDP", "L-FUC", "C16H25N5O15P2", ["S", "R", "R"]),
    ("dTDP-RHA", "dTDP-Rhamnose", ["dTDP-L-rhamnose"],
     "dTDP", "L-RHA", "C16H26N2O15P2", ["R", "R", "S"]),
]


def generate_ndp_sugars(compounds: list[dict]) -> list[dict]:
    """Generate NDP-sugar compounds.

    These are curated with known molecular formulas from biochemistry databases.

    Args:
        compounds: list of all current compounds (for parent lookup)

    Returns:
        list of NDP-sugar compound dicts
    """
    ndp_sugars: list[dict] = []

    for (compound_id, name, aliases, nucleotide, parent_id,
         formula, stereocenters) in CURATED_NDP_SUGARS:

        compound = {
            "id": compound_id,
            "name": name,
            "aliases": aliases,
            "type": "nucleotide_sugar",
            "carbons": 6,  # sugar moiety is always C6 for these
            "chirality": "D" if parent_id.startswith("D-") else "L",
            "formula": formula,
            "stereocenters": stereocenters,
            "modifications": [{"type": "nucleotide", "position": 1}],
            "parent_monosaccharide": parent_id,
            "commercial": False,
            "cost_usd_per_kg": None,
            "metadata": {
                "nucleotide": nucleotide,
                "parent_sugar": parent_id,
            },
            "chebi_id": None,
            "kegg_id": None,
            "pubchem_id": None,
            "inchi": None,
            "smiles": None,
        }
        ndp_sugars.append(compound)

    return ndp_sugars
