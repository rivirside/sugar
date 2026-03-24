"""Generate reactions for lactones.

Supports:
- Lactonization: sugar acid -> lactone (spontaneous, intramolecular dehydration)
- Hydrolysis: lactone -> sugar acid (lactonase, H2O)
"""


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
            {"type": "rule_generated", "rule": f"lactone_{reaction_type}"}
        ],
        "yield": None,
        "cofactor_burden": cofactor_burden,
        "cost_score": 0.94,
        "cofactors": cofactors or [],
        "metadata": {
            "source": "lactone_rule_generation",
            "carbon_count": carbons,
        },
    }


def generate_lactonizations(lactones: list[dict]) -> list[dict]:
    """Generate lactonization/hydrolysis reactions for lactones.

    Each lactone has a parent acid. Generates:
    - acid -> lactone (lactonization, spontaneous or enzyme-catalyzed)
    - lactone -> acid (hydrolysis, lactonase)
    """
    reactions = []

    for lactone in lactones:
        acid_id = lactone.get("metadata", {}).get("parent_acid")
        if not acid_id:
            continue

        carbons = lactone["carbons"]

        # Forward: acid -> lactone (lactonization)
        fwd_id = f"LACT-C{carbons}-{acid_id}-{lactone['id']}"
        reactions.append(
            _base_reaction(
                fwd_id, acid_id, lactone["id"],
                "lactonization", carbons,
            )
        )

        # Reverse: lactone -> acid (hydrolysis)
        rev_id = f"HYDRO-C{carbons}-{lactone['id']}-{acid_id}"
        reactions.append(
            _base_reaction(
                rev_id, lactone["id"], acid_id,
                "hydrolysis", carbons,
            )
        )

    return reactions
