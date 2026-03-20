"""Generate reactions involving phosphosugars."""

from itertools import combinations
from pipeline.reactions.score import compute_cost_score


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


def generate_phosphorylations(phosphosugars: list[dict]) -> list[dict]:
    """Generate phosphorylation reactions: parent monosaccharide -> phosphosugar.

    One reaction per phosphosugar. Requires ATP (cofactor_burden=1.0).
    """
    reactions = []
    for ps in phosphosugars:
        parent_id = ps["parent_monosaccharide"]
        rxn_id = f"PHOS-C{ps['carbons']}-{ps['id']}"
        rxn = _base_reaction(rxn_id, parent_id, ps["id"], "phosphorylation")
        rxn["cofactor_burden"] = 1.0
        rxn["cost_score"] = compute_cost_score(rxn)
        reactions.append(rxn)
    return reactions


def generate_dephosphorylations(phosphosugars: list[dict]) -> list[dict]:
    """Generate dephosphorylation reactions: phosphosugar -> parent monosaccharide.

    Distinct enzyme class from kinases. No cofactor consumed.
    """
    reactions = []
    for ps in phosphosugars:
        parent_id = ps["parent_monosaccharide"]
        rxn_id = f"DEPHOS-C{ps['carbons']}-{ps['id']}"
        rxn = _base_reaction(rxn_id, ps["id"], parent_id, "dephosphorylation")
        reactions.append(rxn)
    return reactions


def generate_mutases(phosphosugars: list[dict]) -> list[dict]:
    """Generate mutase reactions: phosphate migration between positions on same sugar.

    Only between mono-phosphosugars sharing the same parent and stereocenters.
    Bidirectional (forward + reverse).
    """
    reactions = []

    # Group mono-phosphosugars by (parent_monosaccharide, stereocenters_tuple)
    groups: dict[tuple, list[dict]] = {}
    for ps in phosphosugars:
        if len(ps["modifications"]) != 1:
            continue  # Skip bisphosphates
        key = (ps["parent_monosaccharide"], tuple(ps["stereocenters"]))
        groups.setdefault(key, []).append(ps)

    for _key, group in groups.items():
        for sub, prod in combinations(group, 2):
            sub_pos = sub["metadata"]["phosphate_positions"][0]
            prod_pos = prod["metadata"]["phosphate_positions"][0]
            carbons = sub["carbons"]

            # Forward
            fwd_id = f"MUT-C{carbons}-{sub['id']}-{prod_pos}P"
            fwd = _base_reaction(fwd_id, sub["id"], prod["id"], "mutase")
            reactions.append(fwd)

            # Reverse
            rev_id = f"MUT-C{carbons}-{prod['id']}-{sub_pos}P"
            rev = _base_reaction(rev_id, prod["id"], sub["id"], "mutase")
            reactions.append(rev)

    return reactions


def generate_phospho_epimerizations(phosphosugars: list[dict]) -> list[dict]:
    """Generate epimerization reactions between phosphosugars with matching modifications.

    Same logic as standard epimerization: compounds must differ at exactly one
    stereocenter. Additionally, phosphate positions must match exactly.
    """
    reactions = []

    # Group by (parent_type, carbons, phosphate_positions_tuple)
    groups: dict[tuple, list[dict]] = {}
    for ps in phosphosugars:
        key = (
            ps["metadata"]["parent_type"],
            ps["carbons"],
            tuple(sorted(ps["metadata"]["phosphate_positions"])),
        )
        groups.setdefault(key, []).append(ps)

    for _key, group in groups.items():
        for sub, prod in combinations(group, 2):
            if len(sub["stereocenters"]) != len(prod["stereocenters"]):
                continue
            diffs = sum(
                1 for a, b in zip(sub["stereocenters"], prod["stereocenters"]) if a != b
            )
            if diffs != 1:
                continue

            carbons = sub["carbons"]

            # Forward
            fwd_id = f"EPI-C{carbons}-{sub['id']}-{prod['id']}"
            fwd = _base_reaction(fwd_id, sub["id"], prod["id"], "epimerization")
            reactions.append(fwd)

            # Reverse
            rev_id = f"EPI-C{carbons}-{prod['id']}-{sub['id']}"
            rev = _base_reaction(rev_id, prod["id"], sub["id"], "epimerization")
            reactions.append(rev)

    return reactions


def generate_phospho_isomerizations(phosphosugars: list[dict]) -> list[dict]:
    """Generate isomerization reactions between aldose-P and ketose-P.

    Same stereocenter-dropping logic as standard isomerization:
    ketose_stereocenters = aldose_stereocenters[1:]

    Only between phosphosugars with matching phosphate positions.
    """
    reactions = []

    # Separate aldose-derived and ketose-derived phosphosugars
    aldose_ps = [ps for ps in phosphosugars if ps["metadata"]["parent_type"] == "aldose"]
    ketose_ps = [ps for ps in phosphosugars if ps["metadata"]["parent_type"] == "ketose"]

    # Build ketose lookup: (carbons, stereocenters_tuple, phosphate_positions_tuple) -> compound
    ketose_map: dict[tuple, dict] = {}
    for kps in ketose_ps:
        key = (
            kps["carbons"],
            tuple(kps["stereocenters"]),
            tuple(sorted(kps["metadata"]["phosphate_positions"])),
        )
        ketose_map[key] = kps

    for aps in aldose_ps:
        if not aps["stereocenters"]:
            continue  # Skip achiral (C2)

        # Drop first stereocenter to get ketose config
        ketose_config = tuple(aps["stereocenters"][1:])
        phosphate_key = tuple(sorted(aps["metadata"]["phosphate_positions"]))
        lookup = (aps["carbons"], ketose_config, phosphate_key)
        kps = ketose_map.get(lookup)

        if kps is None:
            continue

        carbons = aps["carbons"]

        # Forward: aldose-P -> ketose-P
        fwd_id = f"ISO-C{carbons}-{aps['id']}-{kps['id']}"
        fwd = _base_reaction(fwd_id, aps["id"], kps["id"], "isomerization")
        reactions.append(fwd)

        # Reverse: ketose-P -> aldose-P
        rev_id = fwd_id + "-REV"
        rev = _base_reaction(rev_id, kps["id"], aps["id"], "isomerization")
        reactions.append(rev)

    return reactions
