"""Generate reactions for NDP-sugars.

Supports:
- NDP-sugar epimerization (e.g., UDP-Glucose <-> UDP-Galactose via UDP-glucose 4-epimerase)
"""

from itertools import combinations


def _base_reaction(
    reaction_id: str,
    substrate_id: str,
    product_id: str,
    reaction_type: str,
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
            {"type": "rule_generated", "rule": f"ndp_{reaction_type}"}
        ],
        "yield": None,
        "cofactor_burden": cofactor_burden,
        "cost_score": 0.94,
        "cofactors": cofactors or [],
        "metadata": {
            "source": "ndp_rule_generation",
        },
    }


def generate_ndp_epimerizations(ndp_sugars: list[dict]) -> list[dict]:
    """Generate epimerization reactions between NDP-sugars.

    Two NDP-sugars can epimerize if they have the same nucleotide carrier
    and differ at exactly one stereocenter in the sugar moiety.
    """
    reactions = []

    # Group by nucleotide carrier
    groups: dict[str, list[dict]] = {}
    for c in ndp_sugars:
        nuc = c.get("metadata", {}).get("nucleotide", "")
        if nuc not in groups:
            groups[nuc] = []
        groups[nuc].append(c)

    for nucleotide, members in groups.items():
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

            fwd_id = f"EPI-NDP-{nucleotide}-{a['id']}-{b['id']}"
            rev_id = f"EPI-NDP-{nucleotide}-{b['id']}-{a['id']}"

            reactions.append(
                _base_reaction(fwd_id, a["id"], b["id"], "epimerization")
            )
            reactions.append(
                _base_reaction(rev_id, b["id"], a["id"], "epimerization")
            )

    return reactions
