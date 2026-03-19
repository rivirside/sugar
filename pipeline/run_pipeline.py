"""Pipeline orchestrator: enumerate, validate, generate reactions, and write output."""

import argparse
import json
import os
import sys
from datetime import datetime, timezone

from pipeline.enumerate.monosaccharides import enumerate_all_monosaccharides
from pipeline.enumerate.polyols import generate_polyols
from pipeline.reactions.generate import (
    generate_epimerizations,
    generate_isomerizations,
    generate_reductions,
)
from pipeline.validate.completeness import check_completeness
from pipeline.validate.duplicates import check_duplicates
from pipeline.validate.mass_balance import check_mass_balance

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output")


def _abort(message: str) -> None:
    print(f"[ABORT] {message}", file=sys.stderr)
    sys.exit(1)


def run_pipeline(skip_import: bool = False, refresh: set[str] | None = None) -> dict:
    """Execute the full SUGAR v2 Ring 1 pipeline.

    Args:
        skip_import: If True, skip Ring 2 database import steps.
        refresh: Set of source names to force-refresh cached data for
            (e.g. {"chebi", "kegg", "rhea", "brenda"}). None means use cache.

    Steps:
    1. Enumerate monosaccharides (94 compounds, C2-C7)
    2. Generate polyols with degeneracy detection
    3. Combine all compounds
    4. Validate: completeness and duplicates
    5. Generate reactions: epimerizations, isomerizations, reductions
    6. Mass balance check (ABORT on failure)
    7. Verify reaction ID uniqueness (ABORT on duplicates)
    8. Write output JSON files

    Returns a summary dict with counts and file paths.
    """
    print("=== SUGAR v2 Ring 1 Pipeline ===")

    # Step 1: Enumerate monosaccharides
    print("\n[1/7] Enumerating monosaccharides...")
    monosaccharides = enumerate_all_monosaccharides()
    print(f"  -> {len(monosaccharides)} monosaccharides (aldoses + ketoses, C2-C7)")

    # Step 2: Generate polyols
    print("\n[2/7] Generating polyols...")
    polyols = generate_polyols(monosaccharides)
    print(f"  -> {len(polyols)} polyols (with degeneracy detection)")

    # Step 3: Combine all compounds
    print("\n[3/7] Combining compound sets...")
    all_compounds = monosaccharides + polyols
    print(f"  -> {len(all_compounds)} total compounds")

    # Step 4: Validate
    print("\n[4/7] Validating compound set...")
    completeness_warnings = check_completeness(monosaccharides)
    if completeness_warnings:
        for w in completeness_warnings:
            print(f"  [WARNING] {w}")
    else:
        print("  -> Completeness check passed")

    duplicates = check_duplicates(all_compounds)
    if duplicates:
        for d in duplicates:
            print(f"  [WARNING] Duplicate: {d['original']['id']} / {d['duplicate']['id']}")
    else:
        print("  -> Duplicate check passed")

    # Step 5: Generate reactions
    print("\n[5/7] Generating reactions...")
    epimerizations = generate_epimerizations(all_compounds)
    isomerizations = generate_isomerizations(all_compounds)
    reductions = generate_reductions(all_compounds, polyols)
    all_reactions = epimerizations + isomerizations + reductions
    print(f"  -> {len(epimerizations)} epimerizations")
    print(f"  -> {len(isomerizations)} isomerizations")
    print(f"  -> {len(reductions)} reductions")
    print(f"  -> {len(all_reactions)} total reactions")

    # Step 6: Mass balance check (ABORT on failure)
    print("\n[6/7] Checking mass balance...")
    compound_map = {c["id"]: c for c in all_compounds}
    mass_errors = check_mass_balance(all_reactions, compound_map)
    if mass_errors:
        for e in mass_errors:
            print(f"  [ERROR] {e}", file=sys.stderr)
        _abort(f"Mass balance check failed with {len(mass_errors)} error(s)")
    print("  -> Mass balance check passed")

    # Step 7: Verify reaction ID uniqueness (ABORT on duplicates)
    print("\n[7/7] Verifying reaction ID uniqueness...")
    reaction_ids = [r["id"] for r in all_reactions]
    seen_ids: set[str] = set()
    duplicate_ids = []
    for rid in reaction_ids:
        if rid in seen_ids:
            duplicate_ids.append(rid)
        seen_ids.add(rid)
    if duplicate_ids:
        for did in duplicate_ids:
            print(f"  [ERROR] Duplicate reaction ID: {did}", file=sys.stderr)
        _abort(f"Reaction ID uniqueness check failed: {len(duplicate_ids)} duplicate(s)")
    print("  -> Reaction ID uniqueness check passed")

    import_stats = None
    if not skip_import:
        print("\n=== Ring 2: Database Enrichment ===")
        cache_dir = os.path.join(os.path.dirname(__file__), "cache")
        overrides_path = os.path.join(os.path.dirname(__file__), "data", "match_overrides.json")

        from pipeline.import_.chebi import fetch_chebi_bulk
        from pipeline.import_.kegg import fetch_kegg_compounds_batch
        from pipeline.import_.rhea import fetch_rhea_reactions, classify_reaction_participants
        from pipeline.import_.brenda import fetch_brenda_kinetics
        from pipeline.import_.match import match_all_compounds, load_overrides
        from pipeline.import_.merge import (
            enrich_compound, create_rhea_reaction, find_overlapping_reaction,
            enrich_reaction_with_rhea,
        )
        from pipeline.import_.infer import infer_mirrored_reactions
        from pipeline.validate.mass_balance import check_formula_balance

        refresh_sources = refresh or set()

        # Step R1: Fetch ChEBI data
        print("\n[R1] Fetching ChEBI data...")
        chebi_index = fetch_chebi_bulk(cache_dir, refresh="chebi" in refresh_sources)
        print(f"  -> ChEBI index: {len(chebi_index)} entries")

        # Step R2: Match compounds
        print("\n[R2] Matching compounds to ChEBI...")
        overrides = load_overrides(overrides_path)
        match_report = match_all_compounds(all_compounds, chebi_index, overrides)
        matched = sum(1 for m in match_report.values() if m["chebi_id"])
        print(f"  -> {matched}/{len(all_compounds)} compounds matched")

        # Write match report to cache
        os.makedirs(cache_dir, exist_ok=True)
        match_report_path = os.path.join(cache_dir, "match_report.json")
        with open(match_report_path, "w") as f:
            json.dump(match_report, f, indent=2)

        # Step R3: Enrich compounds
        print("\n[R3] Enriching compounds with external IDs...")
        enriched_compounds = []
        for compound in all_compounds:
            match = match_report.get(compound["id"])
            if match and match["chebi_id"]:
                enriched_compounds.append(enrich_compound(compound, match))
            else:
                enriched_compounds.append(compound)
        all_compounds = enriched_compounds

        # Step R4: Fetch KEGG data for matched compounds
        print("\n[R4] Fetching KEGG data...")
        kegg_ids = [m["kegg_id"] for m in match_report.values() if m.get("kegg_id")]
        if kegg_ids:
            kegg_data = fetch_kegg_compounds_batch(kegg_ids, cache_dir, refresh="kegg" in refresh_sources)
            print(f"  -> {len(kegg_data)} KEGG entries fetched")
        else:
            print("  -> No KEGG IDs to fetch")

        # Step R5: Fetch RHEA reactions
        print("\n[R5] Fetching RHEA reactions...")
        chebi_ids = [m["chebi_id"] for m in match_report.values() if m["chebi_id"]]
        rhea_reactions = fetch_rhea_reactions(chebi_ids, cache_dir, refresh="rhea" in refresh_sources)
        print(f"  -> {len(rhea_reactions)} RHEA reactions found")

        # Step R6: Process RHEA reactions
        print("\n[R6] Processing RHEA reactions...")
        chebi_to_compound = {}
        for compound in all_compounds:
            if compound.get("chebi_id"):
                chebi_to_compound[compound["chebi_id"]] = compound["id"]

        new_reactions = []
        enriched_existing = 0
        for rhea_rxn in rhea_reactions:
            subs = [chebi_to_compound.get(cid) for cid in rhea_rxn["substrate_chebi_ids"] if cid in chebi_to_compound]
            prods = [chebi_to_compound.get(cid) for cid in rhea_rxn["product_chebi_ids"] if cid in chebi_to_compound]
            subs = [s for s in subs if s]
            prods = [p for p in prods if p]

            overlap = find_overlapping_reaction(subs, prods, all_reactions)
            if overlap:
                idx = all_reactions.index(overlap)
                all_reactions[idx] = enrich_reaction_with_rhea(overlap, rhea_rxn)
                enriched_existing += 1
            else:
                new_rxn = create_rhea_reaction(rhea_rxn, chebi_to_compound)
                if new_rxn:
                    new_reactions.append(new_rxn)

        all_reactions.extend(new_reactions)
        print(f"  -> {enriched_existing} existing reactions enriched")
        print(f"  -> {len(new_reactions)} new reactions from RHEA")

        # Step R7: Fetch BRENDA kinetics
        print("\n[R7] Fetching BRENDA kinetics...")
        ec_numbers = list({r.get("ec_number") for r in all_reactions if r.get("ec_number")})
        if ec_numbers:
            brenda_data = fetch_brenda_kinetics(ec_numbers, cache_dir, refresh="brenda" in refresh_sources)
            print(f"  -> {len(brenda_data)} EC numbers with kinetics data")
        else:
            print("  -> No EC numbers to fetch kinetics for")

        # Step R8: Infer D-to-L mirrored reactions (RHEA-sourced only)
        print("\n[R8] Inferring D-to-L mirrored reactions...")
        existing_ids = {r["id"] for r in all_reactions}
        rhea_sourced = [r for r in all_reactions if r.get("rhea_id") or r.get("metadata", {}).get("source") in ("rhea_import", "rhea_discovery")]
        inferred = infer_mirrored_reactions(rhea_sourced, all_compounds, existing_ids)
        all_reactions.extend(inferred)
        print(f"  -> {len(inferred)} inferred mirrored reactions")

        # Step R9: Formula balance check on imported reactions
        print("\n[R9] Checking formula balance on imported reactions...")
        imported_rxns = [r for r in all_reactions if r.get("rhea_id") or r.get("metadata", {}).get("source") == "rhea_import"]
        if imported_rxns:
            compound_map_full = {c["id"]: c for c in all_compounds}
            formula_warnings = check_formula_balance(imported_rxns, compound_map_full)
            for w in formula_warnings:
                print(f"  [WARNING] {w}")
            print(f"  -> {len(imported_rxns)} imported reactions checked, {len(formula_warnings)} warnings")
        else:
            print("  -> No imported reactions to check")

        import_stats = {
            "chebi_matched": matched,
            "rhea_reactions": len(new_reactions),
            "enriched_reactions": enriched_existing,
            "inferred_reactions": len(inferred),
        }

        print("\n=== Ring 2 complete ===")

    # Write output files
    print("\n=== Writing output files ===")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    compounds_path = os.path.join(OUTPUT_DIR, "compounds.json")
    reactions_path = os.path.join(OUTPUT_DIR, "reactions.json")
    metadata_path = os.path.join(OUTPUT_DIR, "pipeline_metadata.json")

    with open(compounds_path, "w") as f:
        json.dump(all_compounds, f, indent=2)
    print(f"  -> {compounds_path}")

    with open(reactions_path, "w") as f:
        json.dump(all_reactions, f, indent=2)
    print(f"  -> {reactions_path}")

    metadata = {
        "pipeline_version": "2.0.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "counts": {
            "monosaccharides": len(monosaccharides),
            "polyols": len(polyols),
            "total_compounds": len(all_compounds),
            "epimerizations": len(epimerizations),
            "isomerizations": len(isomerizations),
            "reductions": len(reductions),
            "total_reactions": len(all_reactions),
        },
        "import_stats": import_stats,
        "completeness_warnings": completeness_warnings,
        "duplicate_warnings": len(duplicates),
        "output_files": {
            "compounds": compounds_path,
            "reactions": reactions_path,
        },
    }

    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2)
    print(f"  -> {metadata_path}")

    # Copy to web/data for Next.js build
    web_data_dir = os.path.join(os.path.dirname(__file__), "..", "web", "data")
    if os.path.exists(web_data_dir):
        import shutil
        shutil.copy2(compounds_path, os.path.join(web_data_dir, "compounds.json"))
        shutil.copy2(reactions_path, os.path.join(web_data_dir, "reactions.json"))
        shutil.copy2(metadata_path, os.path.join(web_data_dir, "pipeline_metadata.json"))
        print(f"  -> Copied to {web_data_dir}")

    print("\n=== Pipeline complete ===")
    print(f"  {len(all_compounds)} compounds, {len(all_reactions)} reactions")

    return metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SUGAR pipeline")
    parser.add_argument("--skip-import", action="store_true", help="Skip database import (Ring 1 only)")
    parser.add_argument("--refresh", action="store_true", help="Force refresh all cached data")
    parser.add_argument("--refresh-chebi", action="store_true", help="Refresh ChEBI cache")
    parser.add_argument("--refresh-kegg", action="store_true", help="Refresh KEGG cache")
    parser.add_argument("--refresh-rhea", action="store_true", help="Refresh RHEA cache")
    parser.add_argument("--refresh-brenda", action="store_true", help="Refresh BRENDA cache")
    args = parser.parse_args()

    refresh_sources = set()
    if args.refresh:
        refresh_sources = {"chebi", "kegg", "rhea", "brenda"}
    else:
        if args.refresh_chebi:
            refresh_sources.add("chebi")
        if args.refresh_kegg:
            refresh_sources.add("kegg")
        if args.refresh_rhea:
            refresh_sources.add("rhea")
        if args.refresh_brenda:
            refresh_sources.add("brenda")

    run_pipeline(skip_import=args.skip_import, refresh=refresh_sources or None)
