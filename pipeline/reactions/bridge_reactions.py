"""Generate bridge reactions connecting derivative compound classes to the main graph.

Creates reactions that link isolated compound types (amino sugars, deoxy sugars,
NDP-sugars) to their parent monosaccharides or phosphosugars, enabling the
pathfinder to route through these intermediates.

Bridge types:
- Amination: monosaccharide -> amino sugar (aminotransferase, glutamine cofactor)
- Deoxygenation: monosaccharide -> deoxy sugar (reductase, NADPH cofactor)
- NDP activation: sugar-1-phosphate -> NDP-sugar (pyrophosphorylase, NTP cofactor)
"""


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
            {"type": "rule_generated", "rule": f"bridge_{reaction_type}"}
        ],
        "yield": None,
        "cofactor_burden": cofactor_burden,
        "cost_score": 0.94,
        "cofactors": cofactors or [],
        "metadata": {
            "source": "bridge_rule_generation",
        },
    }


def generate_amination_bridges(amino_sugars: list[dict]) -> list[dict]:
    """Generate amination/deamination bridges between monosaccharides and amino sugars.

    Each amino sugar with a parent monosaccharide gets:
    - Forward: parent -> amino sugar (aminotransferase, requires glutamine)
    - Reverse: amino sugar -> parent (deaminase)
    """
    reactions = []

    for amino in amino_sugars:
        parent_id = amino.get("parent_monosaccharide")
        if not parent_id:
            continue

        # Only bridge free amino sugars (not N-acetyl, those are reached via nacetylation)
        mods = amino.get("modifications", [])
        if any(m["type"] == "nacetyl" for m in mods):
            continue

        carbons = amino["carbons"]

        # Forward: parent -> amino sugar
        fwd_id = f"AMIN-C{carbons}-{parent_id}-{amino['id']}"
        reactions.append(
            _base_reaction(
                fwd_id, parent_id, amino["id"],
                "transamination", cofactor_burden=1.0,
                cofactors=["glutamine"],
            )
        )

        # Reverse: amino sugar -> parent
        rev_id = f"DEAMIN-C{carbons}-{amino['id']}-{parent_id}"
        reactions.append(
            _base_reaction(
                rev_id, amino["id"], parent_id,
                "hydrolysis",
            )
        )

    return reactions


def generate_deoxygenation_bridges(deoxy_sugars: list[dict]) -> list[dict]:
    """Generate deoxygenation/oxygenation bridges between monosaccharides and deoxy sugars.

    Each deoxy sugar with a parent monosaccharide gets:
    - Forward: parent -> deoxy sugar (reductase, requires NADPH)
    - Reverse: deoxy sugar -> parent (oxygenase)
    """
    reactions = []

    for deoxy in deoxy_sugars:
        parent_id = deoxy.get("parent_monosaccharide")
        if not parent_id:
            continue

        carbons = deoxy["carbons"]

        # Forward: parent -> deoxy sugar
        fwd_id = f"DEOX-C{carbons}-{parent_id}-{deoxy['id']}"
        reactions.append(
            _base_reaction(
                fwd_id, parent_id, deoxy["id"],
                "reduction", cofactor_burden=1.0,
                cofactors=["NADPH"],
            )
        )

        # Reverse: deoxy sugar -> parent
        rev_id = f"REOX-C{carbons}-{deoxy['id']}-{parent_id}"
        reactions.append(
            _base_reaction(
                rev_id, deoxy["id"], parent_id,
                "oxidation", cofactor_burden=1.0,
                cofactors=["NAD+"],
            )
        )

    return reactions


# Map nucleotide type to the NTP cofactor used in pyrophosphorylase reaction
_NTP_COFACTORS = {
    "UDP": "UTP",
    "GDP": "GTP",
    "CMP": "CTP",
    "dTDP": "dTTP",
}


def generate_ndp_activation_bridges(
    ndp_sugars: list[dict],
    all_compounds: list[dict],
) -> list[dict]:
    """Generate NDP-sugar activation bridges.

    Connects sugar-1-phosphates to NDP-sugars via pyrophosphorylase reactions:
    - Forward: sugar-1-P -> NDP-sugar (pyrophosphorylase, requires NTP)
    - Reverse: NDP-sugar -> sugar-1-P (pyrophosphorolysis)

    For NDP-sugars whose parent is not a monosaccharide (e.g., UDP-GlcNAc
    where parent is D-GlcNAc), we look for a sugar-1-phosphate of the parent
    or connect directly to the parent compound.
    """
    reactions = []
    compound_map = {c["id"]: c for c in all_compounds}

    for ndp in ndp_sugars:
        parent_id = ndp.get("parent_monosaccharide")
        nucleotide = ndp.get("metadata", {}).get("nucleotide", "")
        ntp = _NTP_COFACTORS.get(nucleotide, f"{nucleotide}TP")

        if not parent_id:
            continue

        # Try to find the sugar-1-phosphate of the parent
        # Convention: parent-1P (e.g., D-GLC-1P for D-GLC)
        sugar_1p_id = f"{parent_id}-1P"
        sugar_1p = compound_map.get(sugar_1p_id)

        if sugar_1p:
            # Bridge: sugar-1-P <-> NDP-sugar
            fwd_id = f"NDPACT-{nucleotide}-{sugar_1p_id}-{ndp['id']}"
            reactions.append(
                _base_reaction(
                    fwd_id, sugar_1p_id, ndp["id"],
                    "nucleotidyltransfer", cofactor_burden=1.0,
                    cofactors=[ntp],
                )
            )

            rev_id = f"NDPREL-{nucleotide}-{ndp['id']}-{sugar_1p_id}"
            reactions.append(
                _base_reaction(
                    rev_id, ndp["id"], sugar_1p_id,
                    "hydrolysis",
                )
            )
        else:
            # No sugar-1-P exists, bridge directly to parent
            # (simplified, higher cofactor cost)
            parent = compound_map.get(parent_id)
            if parent:
                fwd_id = f"NDPACT-{nucleotide}-{parent_id}-{ndp['id']}"
                reactions.append(
                    _base_reaction(
                        fwd_id, parent_id, ndp["id"],
                        "nucleotidyltransfer", cofactor_burden=2.0,
                        cofactors=[ntp, "ATP"],
                    )
                )

                rev_id = f"NDPREL-{nucleotide}-{ndp['id']}-{parent_id}"
                reactions.append(
                    _base_reaction(
                        rev_id, ndp["id"], parent_id,
                        "hydrolysis",
                    )
                )

    return reactions
