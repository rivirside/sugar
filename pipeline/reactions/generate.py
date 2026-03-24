"""Generate reactions between monosaccharides."""

from itertools import combinations
from pipeline.reactions.score import compute_cost_score


def _reaction_id(reaction_type: str, carbons: int, sub_id: str, prod_id: str) -> str:
    """Generate a descriptive reaction ID like 'EPI-C6-D-GLC-D-MAN'."""
    prefix = {
        "epimerization": "EPI",
        "isomerization": "ISO",
        "reduction": "RED",
    }.get(reaction_type, reaction_type.upper()[:3])
    return f"{prefix}-C{carbons}-{sub_id}-{prod_id}"


def _base_reaction(reaction_id: str, substrate_id: str, product_id: str, reaction_type: str) -> dict:
    """Create a reaction dict with all required fields."""
    rxn = {
        "id": reaction_id,
        "reaction_type": reaction_type,
        "substrates": [substrate_id],
        "products": [product_id],
        "evidence_tier": "hypothetical",
        "evidence_criteria": [],
        "yield": None,
        "cofactor_burden": 0.0,
    }
    rxn["cost_score"] = compute_cost_score(rxn)
    return rxn


def generate_epimerizations(compounds: list[dict]) -> list[dict]:
    """Generate epimerization reactions between compounds differing at exactly one stereocenter.

    Epimerizations are generated within same type (aldose/aldose or ketose/ketose)
    and same carbon count. Both forward and reverse reactions are created.
    """
    reactions = []

    # Group compounds by (type, carbons) — skip derivatives (handled by own modules)
    groups: dict[tuple, list] = {}
    for c in compounds:
        if c["type"] in ("phosphate", "deoxy_sugar", "amino_sugar", "acid", "lactone"):
            continue
        key = (c["type"], c["carbons"])
        groups.setdefault(key, []).append(c)

    for (sugar_type, carbons), group in groups.items():
        for sub, prod in combinations(group, 2):
            # Must have same number of stereocenters
            if len(sub["stereocenters"]) != len(prod["stereocenters"]):
                continue
            # Must differ at exactly one position
            diffs = sum(
                1 for a, b in zip(sub["stereocenters"], prod["stereocenters"]) if a != b
            )
            if diffs != 1:
                continue

            # Forward reaction
            fwd_id = _reaction_id("epimerization", carbons, sub["id"], prod["id"])
            fwd = _base_reaction(fwd_id, sub["id"], prod["id"], "epimerization")
            reactions.append(fwd)

            # Reverse reaction
            rev_id = _reaction_id("epimerization", carbons, prod["id"], sub["id"])
            rev = _base_reaction(rev_id, prod["id"], sub["id"], "epimerization")
            reactions.append(rev)

    return reactions


def generate_isomerizations(compounds: list[dict]) -> list[dict]:
    """Generate isomerization reactions between aldoses and ketoses of the same carbon count.

    An aldose (C-n, stereocenters [c2, c3, ..., cn-1]) isomerizes to the ketose
    whose stereocenters equal the aldose stereocenters with C2 (index 0) dropped:
    ketose_stereocenters = aldose_stereocenters[1:]

    C2 aldoses have no stereocenters and no corresponding ketose, so they are skipped.
    Both forward (aldose->ketose) and reverse reactions are created.
    """
    reactions = []

    # Build lookup: (type, carbons, tuple(stereocenters)) -> compound
    aldoses = [c for c in compounds if c["type"] == "aldose" and c["carbons"] >= 3]
    ketoses = [c for c in compounds if c["type"] == "ketose"]
    ketose_map: dict[tuple, dict] = {
        (c["carbons"], tuple(c["stereocenters"])): c for c in ketoses
    }

    for aldo in aldoses:
        # Ketose stereocenters = aldose stereocenters with C2 (index 0) dropped
        ketose_config = tuple(aldo["stereocenters"][1:])
        ketose = ketose_map.get((aldo["carbons"], ketose_config))
        if ketose is None:
            continue

        # Forward: aldose -> ketose
        fwd_id = _reaction_id("isomerization", aldo["carbons"], aldo["id"], ketose["id"])
        fwd = _base_reaction(fwd_id, aldo["id"], ketose["id"], "isomerization")
        reactions.append(fwd)

        # Reverse: ketose -> aldose
        rev_id = fwd_id + "-REV"
        rev = _base_reaction(rev_id, ketose["id"], aldo["id"], "isomerization")
        reactions.append(rev)

    return reactions


def generate_reductions(compounds: list[dict], polyols: list[dict]) -> list[dict]:
    """Generate reduction reactions from monosaccharides to polyols.

    Each polyol records ALL parent monosaccharides in metadata.reduction_parents.
    A separate reduction reaction is generated for each (parent, polyol) pair.
    Reduction is irreversible and requires NADH (cofactor_burden = 1.0).
    """
    reactions = []

    compound_map = {c["id"]: c for c in compounds}

    for polyol in polyols:
        # Generate one reaction per parent monosaccharide
        all_parents = polyol.get("metadata", {}).get("reduction_parents", [])
        if not all_parents:
            # Fall back to primary parent
            primary = polyol.get("parent_monosaccharide")
            all_parents = [primary] if primary else []

        for parent_id in all_parents:
            parent = compound_map.get(parent_id)
            if parent is None:
                continue

            rxn_id = _reaction_id("reduction", parent["carbons"], parent_id, polyol["id"])
            rxn = _base_reaction(rxn_id, parent_id, polyol["id"], "reduction")
            # Reduction requires NADH
            rxn["cofactor_burden"] = 1.0
            rxn["cost_score"] = compute_cost_score(rxn)
            reactions.append(rxn)

    return reactions
