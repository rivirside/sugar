# Ring 2: Database Enrichment Pipeline

## Overview

Enrich the Ring 1 compound/reaction network (135 compounds, 696 reactions) with data from four external databases: ChEBI, KEGG, RHEA, and BRENDA. This adds real compound identifiers, validated reactions with literature references, enzyme kinetics, and promotes evidence tiers from "hypothetical" to real classifications. New compounds discovered via RHEA (derivatives, nucleotide sugars, etc.) are added to the network immediately.

## Goals

1. Cross-reference enumerated compounds with ChEBI and KEGG IDs
2. Import validated reactions from RHEA (with EC numbers, PMIDs, enzyme names)
3. Add new compounds discovered through RHEA reactions (derivatives, nucleotide sugars, etc.)
4. Pull kinetic parameters (Km, kcat, organism) from BRENDA
5. Promote evidence tiers based on imported data
6. Infer D-to-L mirrored reactions for any validated/predicted reaction found
7. Update frontend to display all new data

## Pipeline Architecture

### Module Structure

```
pipeline/import/
  __init__.py
  chebi.py        # ChEBI: bulk TSV download, fallback to REST per-compound
  kegg.py         # KEGG: async batch REST with caching
  rhea.py         # RHEA: batch SPARQL queries
  brenda.py       # BRENDA: SOAP API per-EC-number
  cache.py        # Shared caching utilities (read/write/refresh/TTL)
  match.py        # Multi-strategy compound matching engine
  merge.py        # Merge imported data into enumerated compounds/reactions
  infer.py        # D-to-L mirroring and gap-fill inference
```

### Pipeline Flow

```
1. Generate (existing)       -> raw compounds + reactions (all hypothetical)
2. Import (new)              -> fetch from ChEBI, KEGG, RHEA, BRENDA
3. Match (new)               -> link enumerated compounds to external IDs
4. Merge (new)               -> enrich compounds/reactions with imported data
5. Add discovered compounds  -> new compounds from RHEA added to network
6. Infer (new)               -> D-to-L mirroring of validated/predicted reactions
7. Update evidence tiers     -> promote reactions based on imported data
8. Re-score (existing)       -> recalculate cost_score with real yields/evidence
9. Validate (existing)       -> mass balance on generated reactions; formula balance on imported reactions
10. Output                   -> enriched compounds.json + reactions.json
```

### CLI

```
python run_pipeline.py                  # Full pipeline (generate + import)
python run_pipeline.py --skip-import    # Generate only (Ring 1 behavior)
python run_pipeline.py --refresh        # Force re-download all cached data
python run_pipeline.py --refresh-chebi  # Refresh specific source
python run_pipeline.py --refresh-kegg
python run_pipeline.py --refresh-rhea
python run_pipeline.py --refresh-brenda
```

## Data Sources

### ChEBI

- **Primary strategy**: Bulk TSV download from ChEBI FTP (compounds table + synonyms table, ~2GB total)
- **Fallback**: If bulk download fails, use ChEBI REST API per-compound (`GET /entity/search/{name}`)
- **Fallback on fallback**: If individual REST call fails, log warning and skip that compound
- **Extracts**: chebi_id, canonical name, synonyms, InChI, SMILES, cross-references (KEGG ID, PubChem ID)
- **Scales**: Bulk download is O(1) regardless of compound count. REST fallback is O(n).

### KEGG

- **Strategy**: Async batch REST requests to `rest.kegg.jp`
  - `/get/{compound_id}` for compound details
  - `/link/reaction/{compound_id}` for linked reactions
- **Rate limiting**: KEGG throttles at ~10 req/sec; batch with delays and exponential backoff
- **Extracts**: kegg_id, pathway memberships, linked reaction IDs
- **Scales**: O(n) with compound count, mitigated by caching

### RHEA

- **Strategy**: Batch SPARQL queries to `sparql.rhea-db.org`
- **Query pattern**: "all reactions where any of these ChEBI IDs appear as substrate or product"
- **Batch size**: ~50 ChEBI IDs per query
- **Extracts**: rhea_id, EC number, substrates/products (as ChEBI IDs), directionality, cross-references to PMIDs
- **New reactions**: RHEA will return reactions not in Ring 1 (enzymatic oxidations, transferases, etc.). These are added as new reaction entries.
- **New compounds**: Reactions involving compounds not in our set cause those compounds to be added (see "Discovered Compounds" section below).

### BRENDA

- **Strategy**: SOAP API per EC number (gathered from RHEA results)
- **Credentials**: Stored in `.env` file (gitignored), loaded at runtime via `python-dotenv`
- **Extracts**: Km (mM), kcat (1/s), specific activity, organism, substrate specificity, temperature/pH optima
- **Scales**: O(e) where e = number of unique EC numbers (~50-200), not compound count

## Compound Matching Engine

### Strategy Order

Applied in sequence; first match wins:

1. **Exact name match** (confidence: high) -- search external DB for the compound's `name` field (e.g., "D-Glucose")
2. **Synonym/alias match** (confidence: high) -- search using entries from `aliases` array (e.g., "Dextrose")
3. **InChI match** (confidence: high) -- generate InChI from our stereochemistry data and match against ChEBI's InChI strings. This is the most reliable structural match. Requires `rdkit` or `pyinchi` for InChI generation from our R/S stereocenters.
4. **Formula match with single-candidate rule** (confidence: medium) -- match on molecular formula. Only applies when formula yields exactly ONE external candidate. If multiple candidates share the formula (e.g., D-Glucose and D-Galactose are both C6H12O6), this strategy produces no match and falls through.
5. **Fuzzy name match** (confidence: low) -- Levenshtein distance on names. Results flagged for manual review, not auto-applied.

### Match Report

Auto-generated at `pipeline/cache/match_report.json`:

```json
{
  "D-GLC": {
    "chebi_id": "CHEBI:17634",
    "kegg_id": "C00031",
    "pubchem_id": "5793",
    "confidence": "high",
    "strategy": "exact_name",
    "chebi_name": "D-glucopyranose"
  },
  "ALDO-C7-RRSRS": {
    "chebi_id": null,
    "kegg_id": null,
    "pubchem_id": null,
    "confidence": null,
    "strategy": "no_match"
  }
}
```

### Match Overrides

A curated `pipeline/data/match_overrides.json` (git-tracked) allows manual corrections:

```json
{
  "D-GLC": {
    "chebi_id": "CHEBI:17634",
    "action": "pin"
  },
  "SOME-COMPOUND": {
    "chebi_id": "CHEBI:99999",
    "action": "reject"
  }
}
```

- `pin`: Force this match regardless of auto-matching result
- `reject`: Explicitly reject a match that auto-matching would accept

The pipeline checks overrides before auto-matching.

## Conflict Resolution

- **Display names**: `pipeline/data/name_mapping.json` remains authoritative. External database names (e.g., ChEBI's "D-glucopyranose") are added to `aliases`, not used as the primary `name`.
- **IDs**: Our compound IDs (e.g., "D-GLC") are never replaced by external IDs. External IDs are stored in dedicated fields (`chebi_id`, `kegg_id`, etc.).
- **Reaction overlap detection**: If RHEA describes a reaction that overlaps with a generated reaction, we enrich the existing reaction rather than creating a duplicate. Overlap is detected by checking whether the RHEA reaction's sugar-class participants (ignoring cofactors) match a generated reaction's substrate/product pair.
- **Reaction IDs for RHEA-sourced reactions**: Use `RHEA:{rhea_id}` as the canonical ID (e.g., `RHEA:12345`). This avoids conflicts with the generated ID scheme (`EPI-C6-...`, `ISO-C6-...`, `RED-C6-...`). When a RHEA reaction enriches an existing generated reaction, the generated ID is kept and the `rhea_id` field is populated.

## Discovered Compounds

When RHEA returns a reaction where one or more participants are not in our compound set:

1. Look up the unknown ChEBI ID to get compound details (name, formula, type)
2. Classify the compound type (phosphate, acid, lactone, nucleotide_sugar, amino_sugar, deoxy_sugar, disaccharide) based on ChEBI ontology classification
3. Add the compound to the network with `parent_monosaccharide` linked if determinable
4. Add the RHEA reaction connecting it to known compounds
5. Tag discovered compounds with `metadata.source: "rhea_discovery"` for traceability

This pulls Ring 3 derivative compounds into the network organically -- only compounds that participate in real reactions with our existing set, not exhaustive enumeration.

### Single-Pass Limitation

RHEA is queried once using the ChEBI IDs of our original compound set. Discovered compounds are NOT re-queried against RHEA for additional reactions. This means reactions between two discovered compounds (where neither was in the original set) will not be captured. This is an acceptable limitation: discovered compounds are included because they connect to our core network, not to exhaustively map their own reaction neighborhoods. Future rings can expand the query frontier.

## Infer Layer: D-to-L Mirroring

For any reaction imported from RHEA with evidence tier "validated" or "predicted":

1. Identify the sugar-class participants in the reaction (compounds with type aldose, ketose, polyol, or derivative types). Cofactors (ATP, NAD+, UDP, etc.) are not sugar-class.
2. For each sugar-class participant, check if a corresponding L/D mirror compound exists in our set.
3. If ALL sugar-class participants have mirrors and the mirrored reaction does not already exist, create it.
4. Cofactors are carried over unchanged to the mirrored reaction (ATP is ATP regardless of D/L context).
5. Tag the inferred reaction with evidence tier "inferred" and evidence_criteria explaining it was mirrored from the validated reaction, referencing the source RHEA ID.

For 1:1 reactions (epimerization, isomerization), this is straightforward -- mirror both compounds. For multi-participant reactions (e.g., `D-Glucose + ATP -> D-Glucose-6-P + ADP`), only the sugar participants are mirrored, cofactors stay the same.

This expands the network with high-confidence inferred reactions.

## Evidence Tier Promotion

After import, reactions are re-classified:

| Condition | Tier |
|-----------|------|
| Found in RHEA with at least one PMID | validated |
| Found in RHEA with EC number but no PMID | predicted |
| Found in RHEA with neither EC number nor PMID (RHEA ID only) | predicted |
| Generated by D-to-L mirroring of a validated/predicted reaction | inferred |
| Rule-based (epimerization/isomerization/reduction) with no database support | hypothetical |

Evidence criteria are recorded as structured data explaining why the tier was assigned:

```json
{
  "evidence_tier": "validated",
  "evidence_criteria": [
    {"source": "rhea", "rhea_id": "RHEA:12345"},
    {"source": "pmid", "ids": ["12345678", "23456789"]},
    {"source": "ec", "ec_number": "5.1.3.3"}
  ]
}
```

## Data Model Changes

### Compound (new fields)

| Field | Type | Source |
|-------|------|--------|
| chebi_id | string or null | ChEBI |
| kegg_id | string or null | ChEBI cross-ref or KEGG |
| pubchem_id | string or null | ChEBI cross-ref |
| inchi | string or null | ChEBI |
| smiles | string or null | ChEBI |

### Reaction (fields already in TypeScript types, now populated)

| Field | Type | Source |
|-------|------|--------|
| ec_number | string or null | RHEA |
| enzyme_name | string or null | RHEA / BRENDA |
| cofactors | string[] | RHEA |
| pmid | string[] | RHEA |
| rhea_id | string or null | RHEA |
| organism | string[] | BRENDA |
| km_mm | number or null | BRENDA |
| kcat_sec | number or null | BRENDA |
| delta_g | number or null | BRENDA (if available) |

## Caching

All cached data lives in `pipeline/cache/` (gitignored).

```
pipeline/cache/
  chebi/
    compounds.tsv           # Bulk download (~2GB)
    names.tsv               # Synonyms table
    index.json              # Parsed/indexed version (much smaller)
  kegg/
    {compound_id}.json      # Per-compound response cache
  rhea/
    query_results.json      # SPARQL query results
  brenda/
    {ec_number}.json        # Per-EC response cache
  match_report.json         # Auto-generated matching results
```

Each cache entry includes a timestamp. The `--refresh` flag bypasses the cache for a full re-download. Source-specific `--refresh-{source}` flags refresh only that source.

## Frontend Updates

### Compound Detail Page
- Add "External IDs" section with clickable links to ChEBI, KEGG, PubChem
- Show InChI string (copyable) and SMILES (copyable)
- Add "Source" badge for discovered compounds (RHEA discovery vs enumerated)

### Reaction Detail Page
- Populate the existing "Ring 2" placeholder section:
  - EC number (linked to BRENDA/ExPASy)
  - Enzyme name
  - Literature references (PMID links to PubMed)
  - RHEA ID (linked to RHEA)
  - Organism list
  - Kinetics: Km, kcat (with units)

### Compound Browser
- Add optional filter: "Has external ID" toggle

### Dashboard
- Add enrichment coverage stats: "X/Y compounds matched to ChEBI", "Z reactions validated"

### TypeScript Type Changes
- Add to Compound interface: `chebi_id`, `kegg_id`, `pubchem_id`, `inchi`, `smiles` (all `string | null`)
- Reaction interface already has all needed optional fields

### Python Pipeline Changes
- The compound dict template in `enumerate/monosaccharides.py` and `enumerate/polyols.py` must be updated to include the new fields (`chebi_id`, `kegg_id`, `pubchem_id`, `inchi`, `smiles`) initialized to `null`. This ensures all compounds have these keys in the JSON output even before enrichment, so the frontend gets `null` rather than `undefined`.
- `run_pipeline.py` must be refactored to accept CLI arguments via `argparse` (`--skip-import`, `--refresh`, etc.). The existing `run_pipeline()` function's return dict is extended with import statistics.

## Mass Balance for Imported Reactions

The existing `check_mass_balance` validates that substrate carbons = product carbons. This works for Ring 1's 1:1 reactions but breaks on RHEA multi-participant reactions where cofactors carry carbons (e.g., `glucose + ATP -> glucose-6-phosphate + ADP` -- ATP and ADP carry carbons the validator doesn't know about).

**Approach**: Two validation modes:
1. **Generated reactions** (epimerization, isomerization, reduction): Continue using the existing carbon-count mass balance. These are always 1:1 and all participants are in our compound map. ABORT on failure.
2. **Imported reactions** (from RHEA): Use molecular formula balance instead. RHEA provides complete participant lists with ChEBI IDs, and ChEBI provides molecular formulas. Validate that the sum of atoms on substrate side equals the sum on product side (accounting for ALL participants including cofactors). Log warnings for mismatches but do not abort -- RHEA is authoritative, so a mismatch likely indicates incomplete formula data rather than an invalid reaction.

The `check_mass_balance` function is extended with a `mode` parameter: `"carbon"` (default, current behavior) or `"formula"` (new, for imported reactions).

## Testing

### Pipeline Tests (pytest)
- Match engine: test each strategy with mock data, test priority ordering, test override application
- ChEBI importer: test bulk parse, test REST fallback, test failure handling
- KEGG importer: test response parsing, test rate limit handling
- RHEA importer: test SPARQL result parsing, test compound discovery
- BRENDA importer: test SOAP response parsing
- Merge: test enrichment of existing compounds, test new compound addition
- Infer: test D-to-L mirroring logic
- Evidence tier promotion: test all tier conditions
- Mass balance: verify imported reactions pass validation

### Frontend Tests (vitest)
- Updated types render correctly
- External ID links generate correct URLs
- New compounds from RHEA display correctly

## Environment

### .env (gitignored)
```
BRENDA_EMAIL=user@example.com
BRENDA_PASSWORD=secret
```

### Dependencies (pipeline)
- `requests` (ChEBI REST, KEGG REST)
- `aiohttp` (async KEGG batching)
- `SPARQLWrapper` (RHEA queries)
- `zeep` (BRENDA SOAP client) -- note: BRENDA has been migrating toward REST; if SOAP endpoint is unavailable, fall back to BRENDA REST API
- `python-dotenv` (credential loading)
- `thefuzz` (fuzzy name matching)
- `rdkit` or `pyinchi` (InChI generation for structural matching)

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| ChEBI bulk download fails or format changes | Automatic fallback to REST API |
| KEGG rate limiting | Exponential backoff + caching |
| RHEA SPARQL endpoint downtime | Cache results; pipeline warns but continues |
| BRENDA credentials invalid | Pipeline skips BRENDA with warning, other sources unaffected |
| False positive compound matches | Multi-strategy confidence scoring + manual override file |
| Imported reactions fail mass balance | Validate all reactions (generated + imported); log failures, don't add invalid reactions |
| Scope creep from discovered compounds | Only add compounds that participate in real reactions with existing compounds |
| BRENDA SOAP endpoint deprecated | Fall back to BRENDA REST API; verify endpoint availability before building against it |
| Mass balance fails on imported reactions | Use formula-level balance (not carbon-count) for RHEA reactions; log warnings, don't abort |
