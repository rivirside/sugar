"""Pipeline orchestrator: enumerate, validate, generate reactions, and write output."""

import argparse
import json
import os
import sys
from datetime import datetime, timezone

from pipeline.enumerate.monosaccharides import enumerate_all_monosaccharides
from pipeline.enumerate.polyols import generate_polyols
from pipeline.enumerate.phosphosugars import generate_phosphosugars
from pipeline.enumerate.deoxy_sugars import generate_deoxy_sugars
from pipeline.enumerate.amino_sugars import generate_amino_sugars
from pipeline.enumerate.sugar_acids import generate_sugar_acids
from pipeline.enumerate.lactones import generate_lactones
from pipeline.enumerate.ndp_sugars import generate_ndp_sugars
from pipeline.reactions.phosphorylation import (
    generate_phosphorylations,
    generate_dephosphorylations,
    generate_mutases,
    generate_phospho_epimerizations,
    generate_phospho_isomerizations,
)
from pipeline.reactions.generate import (
    generate_epimerizations,
    generate_isomerizations,
    generate_reductions,
)
from pipeline.reactions.deoxy_reactions import generate_deoxy_epimerizations
from pipeline.reactions.amino_reactions import (
    generate_amino_epimerizations,
    generate_nacetylations,
)
from pipeline.reactions.acid_reactions import generate_oxidations, generate_acid_epimerizations
from pipeline.reactions.lactone_reactions import generate_lactonizations
from pipeline.reactions.ndp_reactions import generate_ndp_epimerizations
from pipeline.reactions.bridge_reactions import (
    generate_amination_bridges,
    generate_deoxygenation_bridges,
    generate_ndp_activation_bridges,
)
from pipeline.validate.completeness import check_completeness
from pipeline.validate.duplicates import check_duplicates
from pipeline.validate.mass_balance import check_mass_balance
from pipeline.analyze.gap_analysis import run_gap_analysis
from pipeline.analyze.enzyme_index import build_enzyme_index
from pipeline.analyze.tier2_fetch import enrich_enzyme_index

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output")


def _abort(message: str) -> None:
    print(f"[ABORT] {message}", file=sys.stderr)
    sys.exit(1)


def run_pipeline(skip_import: bool = False, refresh: set[str] | None = None) -> dict:
    """Execute the full SUGAR v2 pipeline.

    Args:
        skip_import: If True, skip Ring 2 database import steps.
        refresh: Set of source names to force-refresh cached data for
            (e.g. {"chebi", "kegg", "rhea", "brenda"}). None means use cache.

    Steps:
    1. Enumerate monosaccharides (94 compounds, C2-C7)
    2. Generate polyols with degeneracy detection
    3. Generate phosphosugars (systematic + curated)
    4. Combine all compounds
    5. Validate: completeness and duplicates
    6. Generate reactions (epi, iso, red, phos, dephos, mutase, phospho-epi, phospho-iso)
    7. Mass balance check (ABORT on failure)
    8. Verify reaction ID uniqueness (ABORT on duplicates)

    Returns a summary dict with counts and file paths.
    """
    print("=== SUGAR v2 Pipeline ===")

    # Step 1: Enumerate monosaccharides
    print("\n[1/8] Enumerating monosaccharides...")
    monosaccharides = enumerate_all_monosaccharides()
    print(f"  -> {len(monosaccharides)} monosaccharides (aldoses + ketoses, C2-C7)")

    # Step 2: Generate polyols
    print("\n[2/8] Generating polyols...")
    polyols = generate_polyols(monosaccharides)
    print(f"  -> {len(polyols)} polyols (with degeneracy detection)")

    # Step 3: Generate phosphosugars
    print("\n[3/8] Generating phosphosugars...")
    phosphosugars = generate_phosphosugars(monosaccharides)
    print(f"  -> {len(phosphosugars)} phosphosugars")

    # Step 3b: Generate deoxy sugars
    print("\n[3b/8] Generating deoxy sugars...")
    deoxy_sugars = generate_deoxy_sugars(monosaccharides)
    print(f"  -> {len(deoxy_sugars)} deoxy sugars")

    # Step 3c: Generate amino sugars
    print("\n[3c/8] Generating amino sugars...")
    amino_sugars = generate_amino_sugars(monosaccharides)
    print(f"  -> {len(amino_sugars)} amino sugars")

    # Step 3d: Generate sugar acids
    print("\n[3d/8] Generating sugar acids...")
    sugar_acids = generate_sugar_acids(monosaccharides)
    print(f"  -> {len(sugar_acids)} sugar acids")

    # Step 3e: Generate lactones
    print("\n[3e/8] Generating lactones...")
    lactone_compounds = generate_lactones(sugar_acids)
    print(f"  -> {len(lactone_compounds)} lactones")

    # Step 3f: Generate NDP-sugars
    print("\n[3f/8] Generating NDP-sugars...")
    ndp_sugars = generate_ndp_sugars(monosaccharides + deoxy_sugars + amino_sugars + sugar_acids)
    print(f"  -> {len(ndp_sugars)} NDP-sugars")

    # Step 4: Combine all compounds
    print("\n[4/8] Combining compound sets...")
    all_compounds = (
        monosaccharides + polyols + phosphosugars +
        deoxy_sugars + amino_sugars + sugar_acids + lactone_compounds +
        ndp_sugars
    )
    print(f"  -> {len(all_compounds)} total compounds")

    # Step 5: Validate
    print("\n[5/8] Validating compound set...")
    completeness_warnings = check_completeness(all_compounds)
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

    # Step 6: Generate reactions
    print("\n[6/8] Generating reactions...")
    epimerizations = generate_epimerizations(all_compounds)
    isomerizations = generate_isomerizations(all_compounds)
    reductions = generate_reductions(all_compounds, polyols)
    phosphorylations = generate_phosphorylations(phosphosugars)
    dephosphorylations = generate_dephosphorylations(phosphosugars)
    mutases = generate_mutases(phosphosugars)
    phospho_epimerizations = generate_phospho_epimerizations(phosphosugars)
    phospho_isomerizations = generate_phospho_isomerizations(phosphosugars)
    deoxy_epimerizations = generate_deoxy_epimerizations(deoxy_sugars)
    amino_epimerizations = generate_amino_epimerizations(amino_sugars)
    nacetylations = generate_nacetylations(amino_sugars)
    oxidations = generate_oxidations(sugar_acids)
    acid_epimerizations = generate_acid_epimerizations(sugar_acids)
    lactonizations = generate_lactonizations(lactone_compounds)
    ndp_epimerizations = generate_ndp_epimerizations(ndp_sugars)
    # Bridge reactions connecting compound classes
    amination_bridges = generate_amination_bridges(amino_sugars)
    deoxygenation_bridges = generate_deoxygenation_bridges(deoxy_sugars)
    ndp_activation_bridges = generate_ndp_activation_bridges(ndp_sugars, all_compounds)
    all_reactions = (
        epimerizations + isomerizations + reductions +
        phosphorylations + dephosphorylations + mutases +
        phospho_epimerizations + phospho_isomerizations +
        deoxy_epimerizations +
        amino_epimerizations + nacetylations +
        oxidations + acid_epimerizations + lactonizations +
        ndp_epimerizations +
        amination_bridges + deoxygenation_bridges + ndp_activation_bridges
    )
    print(f"  -> {len(epimerizations)} epimerizations")
    print(f"  -> {len(isomerizations)} isomerizations")
    print(f"  -> {len(reductions)} reductions")
    print(f"  -> {len(phosphorylations)} phosphorylations")
    print(f"  -> {len(dephosphorylations)} dephosphorylations")
    print(f"  -> {len(mutases)} mutases")
    print(f"  -> {len(phospho_epimerizations)} phospho-epimerizations")
    print(f"  -> {len(phospho_isomerizations)} phospho-isomerizations")
    print(f"  -> {len(deoxy_epimerizations)} deoxy-epimerizations")
    print(f"  -> {len(amino_epimerizations)} amino-epimerizations")
    print(f"  -> {len(nacetylations)} N-acetylations/deacetylations")
    print(f"  -> {len(oxidations)} oxidations")
    print(f"  -> {len(acid_epimerizations)} acid-epimerizations")
    print(f"  -> {len(lactonizations)} lactonizations/hydrolyses")
    print(f"  -> {len(ndp_epimerizations)} NDP-sugar epimerizations")
    print(f"  -> {len(amination_bridges)} amination bridges")
    print(f"  -> {len(deoxygenation_bridges)} deoxygenation bridges")
    print(f"  -> {len(ndp_activation_bridges)} NDP-activation bridges")
    print(f"  -> {len(all_reactions)} total reactions")

    # Step 7: Mass balance check (ABORT on failure)
    print("\n[7/8] Checking mass balance...")
    compound_map = {c["id"]: c for c in all_compounds}
    mass_errors = check_mass_balance(all_reactions, compound_map)
    if mass_errors:
        for e in mass_errors:
            print(f"  [ERROR] {e}", file=sys.stderr)
        _abort(f"Mass balance check failed with {len(mass_errors)} error(s)")
    print("  -> Mass balance check passed")

    # Step 8: Verify reaction ID uniqueness (ABORT on duplicates)
    print("\n[8/8] Verifying reaction ID uniqueness...")
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

    # === Ring 4: Enzyme Gap Analysis ===
    print("\n=== Ring 4: Enzyme Gap Analysis ===")

    # Step G1: Build enzyme index
    print("\n[G1] Building enzyme index...")
    enzyme_index = build_enzyme_index(all_reactions)
    print(f"  -> {len(enzyme_index)} EC families indexed")

    # Step G1.5: Tier 2 enrichment (UniProt/PDB)
    if not skip_import:
        print("\n[G1.5] Enriching enzyme index with UniProt/PDB data...")
        enzyme_index = enrich_enzyme_index(enzyme_index, use_cache=True)
        tier2_count = sum(1 for e in enzyme_index.values() if e.get("family_size") is not None)
        print(f"  -> {tier2_count}/{len(enzyme_index)} EC families enriched with Tier 2 data")
    else:
        print("\n[G1.5] Skipping Tier 2 enrichment (--skip-import)")

    # Step G2: Run gap analysis
    print("\n[G2] Running gap analysis...")
    all_reactions, gap_metadata = run_gap_analysis(
        all_compounds, all_reactions, enzyme_index=enzyme_index
    )
    print(f"  -> {gap_metadata['reactions_analyzed']} reactions analyzed")
    print(f"  -> {gap_metadata['coverage_direct']} direct enzyme matches")
    print(f"  -> {gap_metadata['coverage_cross_substrate']} cross-substrate candidates")
    print(f"  -> {gap_metadata['coverage_family_only']} EC family only")
    print(f"  -> {gap_metadata['coverage_none']} no coverage")
    print(f"  -> avg engineerability: {gap_metadata['avg_engineerability_score']:.4f}")

    print("\n=== Ring 4 complete ===")

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

    # Export enzyme index
    enzyme_index_path = os.path.join(OUTPUT_DIR, "enzyme_index.json")
    with open(enzyme_index_path, "w") as f:
        json.dump(enzyme_index, f, indent=2)
    print(f"  -> {enzyme_index_path}")

    metadata = {
        "pipeline_version": "2.0.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "counts": {
            "monosaccharides": len(monosaccharides),
            "polyols": len(polyols),
            "phosphosugars": len(phosphosugars),
            "deoxy_sugars": len(deoxy_sugars),
            "amino_sugars": len(amino_sugars),
            "sugar_acids": len(sugar_acids),
            "lactones": len(lactone_compounds),
            "ndp_sugars": len(ndp_sugars),
            "total_compounds": len(all_compounds),
            "epimerizations": len(epimerizations),
            "isomerizations": len(isomerizations),
            "reductions": len(reductions),
            "phosphorylations": len(phosphorylations),
            "dephosphorylations": len(dephosphorylations),
            "mutases": len(mutases),
            "phospho_epimerizations": len(phospho_epimerizations),
            "phospho_isomerizations": len(phospho_isomerizations),
            "deoxy_epimerizations": len(deoxy_epimerizations),
            "amino_epimerizations": len(amino_epimerizations),
            "nacetylations": len(nacetylations),
            "oxidations": len(oxidations),
            "acid_epimerizations": len(acid_epimerizations),
            "lactonizations": len(lactonizations),
            "ndp_epimerizations": len(ndp_epimerizations),
            "total_reactions": len(all_reactions),
        },
        "gap_analysis": gap_metadata,
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
        shutil.copy2(enzyme_index_path, os.path.join(web_data_dir, "enzyme_index.json"))
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
