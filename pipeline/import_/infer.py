"""D-to-L mirroring inference."""

from pipeline.reactions.score import compute_cost_score

SUGAR_TYPES = {"aldose", "ketose", "polyol", "phosphate", "acid", "lactone", "amino_sugar", "nucleotide_sugar", "deoxy_sugar", "disaccharide"}


def find_mirror_compound(compound_id: str, compounds: list[dict]) -> dict | None:
    compound_map = {c["id"]: c for c in compounds}
    compound = compound_map.get(compound_id)
    if not compound or compound["chirality"] not in ("D", "L"):
        return None
    mirrored_centers = ["S" if s == "R" else "R" for s in compound["stereocenters"]]
    for c in compounds:
        if (c["id"] != compound_id and c["type"] == compound["type"] and c.get("carbons") == compound.get("carbons") and c["stereocenters"] == mirrored_centers):
            return c
    return None


def infer_mirrored_reactions(reactions: list[dict], compounds: list[dict], existing_reaction_ids: set[str]) -> list[dict]:
    compound_map = {c["id"]: c for c in compounds}
    inferred = []
    for rxn in reactions:
        if rxn.get("evidence_tier") not in ("validated", "predicted"):
            continue
        sugar_substrates, cofactor_substrates = [], []
        for sid in rxn["substrates"]:
            comp = compound_map.get(sid)
            if comp and comp["type"] in SUGAR_TYPES:
                sugar_substrates.append(sid)
            else:
                cofactor_substrates.append(sid)
        sugar_products, cofactor_products = [], []
        for pid in rxn["products"]:
            comp = compound_map.get(pid)
            if comp and comp["type"] in SUGAR_TYPES:
                sugar_products.append(pid)
            else:
                cofactor_products.append(pid)

        mirrored_substrates, all_have_mirrors = [], True
        for sid in sugar_substrates:
            mirror = find_mirror_compound(sid, compounds)
            if mirror:
                mirrored_substrates.append(mirror["id"])
            else:
                all_have_mirrors = False
                break
        if not all_have_mirrors:
            continue

        mirrored_products = []
        for pid in sugar_products:
            mirror = find_mirror_compound(pid, compounds)
            if mirror:
                mirrored_products.append(mirror["id"])
            else:
                all_have_mirrors = False
                break
        if not all_have_mirrors:
            continue

        new_substrates = mirrored_substrates + cofactor_substrates
        new_products = mirrored_products + cofactor_products
        source_rhea = rxn.get("rhea_id", rxn["id"])
        new_id = f"INFER-{source_rhea}-L"
        if new_id in existing_reaction_ids:
            continue

        mirror_rxn = {
            "id": new_id, "reaction_type": rxn.get("reaction_type", "isomerization"),
            "substrates": new_substrates, "products": new_products,
            "evidence_tier": "inferred",
            "evidence_criteria": [{"source": "d_to_l_mirror", "source_reaction": source_rhea}],
            "yield": None, "cofactor_burden": rxn.get("cofactor_burden", 0.0),
            "ec_number": rxn.get("ec_number"), "enzyme_name": rxn.get("enzyme_name"),
            "cofactors": rxn.get("cofactors", []), "pmid": [], "rhea_id": None,
            "organism": [], "km_mm": None, "kcat_sec": None, "delta_g": None,
            "metadata": {"source": "d_to_l_inference", "source_reaction": source_rhea},
        }
        mirror_rxn["cost_score"] = compute_cost_score(mirror_rxn)
        inferred.append(mirror_rxn)
    return inferred
