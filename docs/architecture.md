# Architecture Guide

How SUGAR is built and how to extend it.

## Overview

SUGAR has two main components:

1. **Pipeline** (`pipeline/`): A Python program that generates compound and reaction data from stereochemistry rules, optionally enriches it with external databases, validates it, and writes JSON output files.
2. **Frontend** (`web/`): A Next.js application that loads the JSON output at build time and provides client-side pathway finding, compound/reaction browsing, and network visualization.

There is no runtime server. The pipeline runs offline to produce static data, and the frontend operates entirely in the browser.

## The Ring model

The pipeline is organized into concentric rings. Each ring adds a layer of data without modifying the layers beneath it.

```
Ring 4 (planned): disaccharides, advanced scoring
Ring 3 (planned): acids, lactones, amino sugars, nucleotide sugars, deoxy sugars
Ring 2: database enrichment (ChEBI, KEGG, RHEA, BRENDA)
Ring 1: core compounds and reactions from first principles
```

**Ring 1** is the foundation. It enumerates all possible stereoisomers of C2-C7 monosaccharides, generates their polyol reduction products and phosphorylated derivatives, then creates reactions between them using structural comparison rules. Every reaction starts with the "hypothetical" evidence tier because it is generated from theory, not observed in a database.

**Ring 2** enriches Ring 1 data. It fetches external database records, matches them to Ring 1 compounds and reactions, and upgrades evidence tiers when matches are found. A compound matched in ChEBI gets a `chebi_id`. A reaction matched in RHEA gets upgraded from "hypothetical" to "validated." Ring 2 is optional and can be skipped with `--skip-import`.

**Rings 3 and 4** are planned extensions. Ring 3 will add more sugar derivative classes. Ring 4 will add multi-sugar compounds and advanced scoring.

## Pipeline stages

The pipeline runs in numbered steps, printed to the console as it executes.

### Ring 1 steps

**Step 1: Enumerate monosaccharides** (`enumerate/monosaccharides.py`)

Generates all aldose and ketose stereoisomers for carbon lengths 2 through 7. Aldoses have their carbonyl at C1 and produce 2^(n-2) stereoisomers for a sugar with n carbons. Ketoses have their carbonyl at C2 and produce 2^(n-3) stereoisomers. Each compound gets a systematic ID, stereocenter array, molecular formula, and human-readable name from `data/name_mapping.json` when available.

Output: 94 compounds (63 aldoses + 31 ketoses).

**Step 2: Generate polyols** (`enumerate/polyols.py`)

Reduces each monosaccharide to its corresponding polyol (sugar alcohol). The carbonyl group becomes a hydroxyl, which can create a new stereocenter or, in some cases, produce an achiral molecule. The key challenge is degeneracy: multiple monosaccharides can reduce to the same polyol. For example, D-Glucose and L-Gulose both produce Sorbitol.

Degeneracy is detected by computing a canonical configuration key: `max(config_string, reversed_config_string)`. If two monosaccharides produce the same canonical key, they map to the same polyol. The polyol tracks all of its parents in metadata.

Output: 41 unique polyols.

**Step 3: Generate phosphosugars** (`enumerate/phosphosugars.py`)

Creates phosphorylated derivatives in two modes:

- Systematic: all C6 aldohexose and ketohexose stereoisomers combined with standard phosphorylation positions (C1, C3, C4, C6 for mono-phosphates; C1,6 and C3,6 for bis-phosphates on aldohexoses; C1,6 for bis-phosphates on ketohexoses). This produces 136 compounds.
- Curated: 8 biologically important phosphosugars from other carbon lengths (G3P, DHAP, E4P, R5P, Ru5P, Xu5P, S7P, F26BP).

Each phosphate ester adds H3PO4 minus H2O to the parent formula (net +1P, +3O, +1H per phosphate group).

Output: 144 phosphosugars.

**Step 4: Combine and validate**

All compound sets are merged (94 + 41 + 144 = 279 compounds). The completeness validator checks that each expected stereoisomer group has the correct count. The duplicate validator checks for identical stereocenters within the same compound group.

**Step 5: Generate reactions** (`reactions/generate.py`, `reactions/phosphorylation.py`)

Eight reaction types are generated:

| Type | Logic | Count |
|------|-------|-------|
| Epimerization | Two compounds of the same type and carbon count differ at exactly one stereocenter | 478 |
| Isomerization | An aldose and a ketose of the same carbon count where dropping the aldose C2 stereocenter yields the ketose stereocenters | 124 |
| Reduction | A monosaccharide reduces to its polyol (cofactor: NADH) | 94 |
| Phosphorylation | A monosaccharide gains a phosphate group (cofactor: ATP) | 144 |
| Dephosphorylation | A phosphosugar loses its phosphate group(s) | 144 |
| Mutase | A phosphate migrates from one position to another on the same sugar | 288 |
| Phospho-epimerization | Two phosphosugars with the same phosphate positions differ at exactly one stereocenter | 506 |
| Phospho-isomerization | Aldose-P converts to ketose-P (or vice versa) with the same phosphate positions | 162 |

All reactions are bidirectional (A to B and B to A counted separately) except reductions, which are irreversible.

**Step 6: Score and validate reactions** (`reactions/score.py`, `validate/mass_balance.py`)

Each reaction gets a cost score: `W1*(1-yield) + W2*cofactor_burden + W3*evidence_penalty + W4`. Mass balance is checked (substrate carbons must equal product carbons). The pipeline aborts if any reaction fails mass balance.

### Ring 2 steps (optional)

Ring 2 runs 9 additional steps that fetch external data, match it to Ring 1 compounds and reactions, and merge enrichment fields. These steps are skipped when `--skip-import` is passed.

The matching engine (`import_/match.py`) uses five strategies in priority order: override pin, override reject, exact name, alias match, formula unique match, and fuzzy name match. Results are cached locally to avoid repeated API calls.

### Output

The pipeline writes three files to `pipeline/output/`:

- `compounds.json`: array of all compounds
- `reactions.json`: array of all reactions
- `pipeline_metadata.json`: generation timestamp, version, counts, warnings

It then copies these files to `web/data/` so the frontend can import them at build time.

## Directory layout

```
pipeline/
  run_pipeline.py                 Main orchestrator, runs all steps
  enumerate/
    monosaccharides.py            Aldose/ketose stereoisomer generation
    polyols.py                    Polyol generation with degeneracy detection
    phosphosugars.py              Phosphosugar enumeration (systematic + curated)
  reactions/
    generate.py                   Epimerization, isomerization, reduction
    phosphorylation.py            Phospho reactions (phos, dephos, mutase, epi, iso)
    score.py                      Cost scoring formula
  validate/
    completeness.py               Expected stereoisomer counts per group
    duplicates.py                 Identical-stereocenter detection
    mass_balance.py               Carbon/formula balance checks (fatal on failure)
  import_/
    chebi.py                      ChEBI bulk data fetcher
    kegg.py                       KEGG compound fetcher
    rhea.py                       RHEA reaction fetcher (SPARQL)
    brenda.py                     BRENDA kinetics fetcher (SOAP API)
    match.py                      Multi-strategy compound matching
    merge.py                      Enrichment field merger
    infer.py                      D-to-L reaction inference
    cache.py                      Local caching for API responses
  data/
    name_mapping.json             Stereo-config to human-readable name mapping
    match_overrides.json          Manual match pins and rejects
  output/                         Generated JSON files
  cache/                          Cached API responses (gitignored)
  tests/                          pytest test suite (116 tests)

web/
  app/
    page.tsx                      Dashboard with pathway finder entry
    pathways/page.tsx             Pathway finder (Yen's K-shortest paths)
    compounds/page.tsx            Compound browser with filters
    reactions/page.tsx            Reaction browser with filters
    compound/[id]/page.tsx        Compound detail page
    reaction/[id]/page.tsx        Reaction detail page
    network/page.tsx              Cytoscape network graph
    about/page.tsx                About page
    layout.tsx                    Root layout (navbar, command palette)
  components/
    nav-bar.tsx                   Top navigation
    command-palette.tsx           Cmd+K global search
    compound-search.tsx           Typeahead compound search (Fuse.js)
    evidence-filter.tsx           Evidence tier toggle (URL-synced)
    evidence-badge.tsx            Evidence tier badge
    pathway-list.tsx              K-shortest paths list
    pathway-detail.tsx            Step-by-step pathway breakdown
    stat-card.tsx                 Dashboard statistics cards
    ui/                           Reusable UI primitives (shadcn)
  lib/
    types.ts                      TypeScript interfaces
    data.ts                       JSON data loader, compound/reaction maps
    pathfinding.ts                Dijkstra + Yen's K-shortest paths
    graph.ts                      Adjacency list builder
    search.ts                     Fuse.js search index
    evidence-filter.ts            Evidence tier filter state
    utils.ts                      Color maps, formatting helpers
  data/                           Pipeline output (copied during build)
```

## How to add a new compound class

To add a new type of compound (e.g., amino sugars), follow this pattern:

1. **Create an enumeration module** at `pipeline/enumerate/your_type.py`. Export a function `generate_your_type(compounds: list[dict]) -> list[dict]` that takes existing compounds (so you can derive from them) and returns new compound dicts. Each compound must have all required fields: `id`, `name`, `type`, `carbons`, `chirality`, `formula`, `stereocenters`, `modifications`, `parent_monosaccharide`, `metadata`, and the external ID fields.

2. **Write tests first** at `pipeline/tests/test_your_type.py`. Test counts, field correctness, formula accuracy, and edge cases.

3. **Create a reaction module** if your compound class has its own reaction types. Place it at `pipeline/reactions/your_reactions.py`. Each generator function should take `compounds: list[dict]` and return `list[dict]` of reactions.

4. **Update validation** in `pipeline/validate/completeness.py` to add expected counts for your new compound groups. Update `duplicates.py` if your grouping key needs a new discriminator.

5. **Guard existing generators** in `pipeline/reactions/generate.py`. If existing generators (epimerization, isomerization) should not process your new type, add a `continue` filter for `c["type"] == "your_type"`.

6. **Wire into the pipeline** in `pipeline/run_pipeline.py`: import your generator, add an enumeration step, add reaction generation calls, and update the metadata counts.

7. **Update the frontend types** in `web/lib/types.ts` if needed (add to `CompoundType` union, add new `ReactionType` values).

## How to add a new reaction type

1. **Define the rule** clearly: what structural relationship between two compounds qualifies them for this reaction? Write it down before coding.

2. **Implement the generator** function. It should iterate over all relevant compound pairs, check the structural relationship, and emit reaction dicts. Use the existing `_base_reaction` helper pattern (see `reactions/generate.py` or `reactions/phosphorylation.py`).

3. **Choose the right module**: if this reaction type only involves a specific compound class (like phosphosugars), put it in that class's reaction module. If it spans compound classes, put it in `generate.py`.

4. **Test thoroughly**: test that the count is correct, that specific known reactions appear, and that invalid pairs are excluded.

5. **Wire it in** by calling your generator from `run_pipeline.py` and adding its count to `pipeline_metadata.json`.

## Testing strategy

**Pipeline tests** (`pipeline/tests/`, pytest):

- Unit tests for each enumeration module (correct counts, field values, edge cases)
- Unit tests for each reaction generator (correct pairs, exclusions, bidirectionality)
- Scoring tests (formula correctness, evidence tier defaults)
- Validation tests (completeness, duplicates, mass balance)
- Import/enrichment tests (matching strategies, merge logic, inference)
- Integration test (`test_pipeline_integration.py`) that runs the full pipeline and checks aggregate counts

**Frontend tests** (`web/__tests__/`, Vitest):

- Pathfinding algorithm tests (Dijkstra, Yen's K-shortest paths)
- Edge cases (unreachable targets, max steps, cost computation)

All tests are deterministic and require no network access (Ring 2 import tests mock API responses).

## Frontend tech stack

- Next.js 16 with App Router
- React 19
- TypeScript 5
- Tailwind CSS 4 with shadcn UI components
- Cytoscape.js for network graph visualization
- Fuse.js for fuzzy compound search
- Vitest for testing

The frontend is entirely client-side. Data is imported as JSON at build time. Pathfinding runs in the browser. No API calls are made at runtime.
