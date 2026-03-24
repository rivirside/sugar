"""Generate reactions for sugar acids.

Supports:
- Oxidation: monosaccharide -> sugar acid (NAD+/NADP+ dependent)
- Epimerization between sugar acids with same acid type
"""

from itertools import combinations


def _base_reaction(
    reaction_id: str,
    substrate_id: str,
    product_id: str,
    reaction_type: str,
    carbons: int,
    cofactor_burden: float = 0.0,
    cofactors: list[str] | None = None,
) -> dict:
    return {
        "id": reaction_id,
        "substrates": [substrate_id],
        "products": [product_id],
        "reaction_type": reaction_type,
        "evidence_tier": "hypothetical",
        "evidence_criteria": [
            {"type": "rule_generated", "rule": f"acid_{reaction_type}"}
        ],
        "yield": None,
        "cofactor_burden": cofactor_burden,
        "cost_score": 0.94,
        "cofactors": cofactors or [],
        "metadata": {
            "source": "acid_rule_generation",
            "carbon_count": carbons,
        },
    }


def generate_oxidations(sugar_acids: list[dict]) -> list[dict]:
    """Generate oxidation reactions: monosaccharide -> sugar acid.

    Each sugar acid has a parent monosaccharide. Oxidation requires NAD+.
    """
    reactions = []

    for acid in sugar_acids:
        parent_id = acid.get("parent_monosaccharide")
        if not parent_id:
            continue

        carbons = acid["carbons"]
        acid_type = acid.get("metadata", {}).get("acid_type", "")

        # Forward: parent -> acid (oxidation, NAD+ dependent)
        fwd_id = f"OX-{acid_type.upper()}-C{carbons}-{parent_id}-{acid['id']}"
        reactions.append(
            _base_reaction(
                fwd_id, parent_id, acid["id"],
                "oxidation", carbons,
                cofactor_burden=1.0,
                cofactors=["NAD+"],
            )
        )

    return reactions


def generate_acid_epimerizations(sugar_acids: list[dict]) -> list[dict]:
    """Generate epimerization reactions between sugar acids.

    Two sugar acids can epimerize if they have the same carbon count,
    same acid type (both uronic or both aldonic), and differ at exactly
    one stereocenter.
    """
    reactions = []

    # Group by (carbons, acid_type)
    groups: dict[tuple, list[dict]] = {}
    for c in sugar_acids:
        acid_type = c.get("metadata", {}).get("acid_type", "")
        key = (c["carbons"], acid_type)
        if key not in groups:
            groups[key] = []
        groups[key].append(c)

    for (carbons, acid_type), members in groups.items():
        if len(members) < 2:
            continue

        for a, b in combinations(members, 2):
            stereo_a = a["stereocenters"]
            stereo_b = b["stereocenters"]

            if len(stereo_a) != len(stereo_b):
                continue

            diffs = sum(1 for sa, sb in zip(stereo_a, stereo_b) if sa != sb)
            if diffs != 1:
                continue

            fwd_id = f"EPI-ACID-C{carbons}-{a['id']}-{b['id']}"
            rev_id = f"EPI-ACID-C{carbons}-{b['id']}-{a['id']}"

            reactions.append(
                _base_reaction(fwd_id, a["id"], b["id"], "epimerization", carbons)
            )
            reactions.append(
                _base_reaction(rev_id, b["id"], a["id"], "epimerization", carbons)
            )

    return reactions
