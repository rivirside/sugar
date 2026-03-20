# Data Guide

How SUGAR's data is structured, where it comes from, and the biochemistry behind it.

## Compound data model

Every compound in `compounds.json` has these fields:

```json
{
  "id": "D-GLC",
  "name": "D-Glucose",
  "aliases": ["Dextrose"],
  "type": "aldose",
  "carbons": 6,
  "chirality": "D",
  "formula": "C6H12O6",
  "stereocenters": ["R", "S", "S", "R"],
  "modifications": null,
  "parent_monosaccharide": null,
  "commercial": false,
  "cost_usd_per_kg": null,
  "metadata": {},
  "chebi_id": null,
  "kegg_id": null,
  "pubchem_id": null,
  "inchi": null,
  "smiles": null
}
```

### Field reference

| Field | Type | Description |
|-------|------|-------------|
| `id` | string | Unique systematic identifier. Aldoses: `D-GLC`. Ketoses: `D-FRU`. Polyols: `D-SORBITOL`. Phosphosugars: `D-GLC-6P`. |
| `name` | string | Human-readable name. Known sugars use their common name; unknown stereoisomers use systematic names like "D-aldohexose-RSRS". |
| `aliases` | string[] | Alternative names (e.g., "Dextrose" for D-Glucose). Empty array if none. |
| `type` | string | One of: `aldose`, `ketose`, `polyol`, `phosphate`. Future rings will add `acid`, `lactone`, `amino_sugar`, `nucleotide_sugar`, `deoxy_sugar`, `disaccharide`. |
| `carbons` | number | Carbon chain length (2 through 7 for Ring 1 compounds). |
| `chirality` | string | `"D"`, `"L"`, or `"achiral"`. Based on the configuration of the highest-numbered stereocenter. |
| `formula` | string | Molecular formula in Hill notation (e.g., `"C6H12O6"`). |
| `stereocenters` | string[] | Array of R/S designations for each stereocenter, ordered from lowest to highest carbon number. Empty for achiral compounds. |
| `modifications` | array or null | For phosphosugars: `[{"type": "phosphate", "position": 6}]`. Null for unmodified compounds. |
| `parent_monosaccharide` | string or null | ID of the parent compound for polyols and phosphosugars. Null for monosaccharides. |
| `commercial` | boolean | Whether the compound is commercially available (populated by Ring 2). |
| `cost_usd_per_kg` | number or null | Commercial cost if known. |
| `metadata` | object | Additional data. Contents vary by compound type (see below). |
| `chebi_id` | string or null | ChEBI database identifier (populated by Ring 2). |
| `kegg_id` | string or null | KEGG compound identifier (populated by Ring 2). |
| `pubchem_id` | string or null | PubChem compound identifier (populated by Ring 2). |
| `inchi` | string or null | InChI structural identifier (populated by Ring 2). |
| `smiles` | string or null | SMILES notation (populated by Ring 2). |

### Metadata by compound type

**Polyols** include:
- `reduction_parents`: array of parent monosaccharide IDs (multiple when degenerate)
- `degeneracy`: number of distinct parents that reduce to this polyol

**Phosphosugars** include:
- `phosphate_positions`: array of position numbers (e.g., `[1, 6]` for a 1,6-bisphosphate)
- `parent_type`: the type of the parent compound (`"aldose"` or `"ketose"`)
- `curated`: boolean, true for the 8 hand-picked biologically important phosphosugars

## Reaction data model

Every reaction in `reactions.json` has these fields:

```json
{
  "id": "EPI-C6-D-GLC-D-MAN",
  "reaction_type": "epimerization",
  "substrates": ["D-GLC"],
  "products": ["D-MAN"],
  "evidence_tier": "hypothetical",
  "evidence_criteria": [],
  "yield": null,
  "cofactor_burden": 0.0,
  "cost_score": 0.94
}
```

### Field reference

| Field | Type | Description |
|-------|------|-------------|
| `id` | string | Descriptive ID encoding the reaction type and participants. |
| `reaction_type` | string | One of: `epimerization`, `isomerization`, `reduction`, `phosphorylation`, `dephosphorylation`, `mutase`, `phospho_epimerization`, `phospho_isomerization`. |
| `substrates` | string[] | Array of substrate compound IDs (currently always length 1). |
| `products` | string[] | Array of product compound IDs (currently always length 1). |
| `evidence_tier` | string | One of: `validated`, `predicted`, `inferred`, `hypothetical`. See "Evidence tiers" below. |
| `evidence_criteria` | array | Supporting evidence entries. Empty for hypothetical reactions. Populated by Ring 2 with database references. |
| `yield` | number or null | Expected reaction yield (0.0 to 1.0). Null means unknown. |
| `cofactor_burden` | number | Metabolic cost of required cofactors. 0.0 for no cofactor, 1.0 for ATP or NADH. |
| `cost_score` | number | Composite cost used for pathfinding. Lower is better. |

### Reaction ID format

Each reaction type uses a consistent ID pattern:

| Type | Pattern | Example |
|------|---------|---------|
| Epimerization | `EPI-C{n}-{substrate}-{product}` | `EPI-C6-D-GLC-D-MAN` |
| Isomerization | `ISO-C{n}-{aldose}-{ketose}` | `ISO-C6-D-GLC-D-FRU` |
| Reduction | `RED-C{n}-{sugar}-{polyol}` | `RED-C6-D-GLC-D-SORBITOL` |
| Phosphorylation | `PHOS-C{n}-{phosphosugar}` | `PHOS-C6-D-GLC-6P` |
| Dephosphorylation | `DEPHOS-C{n}-{phosphosugar}` | `DEPHOS-C6-D-GLC-6P` |
| Mutase | `MUT-C{n}-{sugar}-{pos1}P-{pos2}P` | `MUT-C6-D-GLC-1P-6P` |
| Phospho-epimerization | `EPI-C{n}-{substrate}-{product}` | `EPI-C6-D-GLC-6P-D-MAN-6P` |
| Phospho-isomerization | `ISO-C{n}-{aldose-P}-{ketose-P}` | `ISO-C6-D-GLC-6P-D-FRU-6P` |

### Cost scoring formula

```
cost_score = W1 * (1 - yield) + W2 * cofactor_burden + W3 * evidence_penalty + W4
```

Where:
- W1 = 1.0 (yield loss weight)
- W2 = 0.5 (cofactor burden weight)
- W3 = 0.3 (evidence uncertainty weight)
- W4 = 0.2 (base step penalty, constant per reaction)

Evidence penalties:

| Tier | Default yield | Penalty |
|------|--------------|---------|
| validated | 1.0 | 0.0 |
| predicted | 0.8 | 0.1 |
| inferred | 0.6 | 0.3 |
| hypothetical | 0.5 | 0.8 |

A hypothetical reaction with no cofactor scores: 1.0*(1-0.5) + 0.5*0.0 + 0.3*0.8 + 0.2 = 0.94.

A hypothetical reduction (cofactor: NADH = 1.0) scores: 1.0*(1-0.5) + 0.5*1.0 + 0.3*0.8 + 0.2 = 1.44.

Lower cost means the pathfinding algorithm prefers that reaction.

## Evidence tiers

Evidence tiers indicate how confident we are that a reaction exists in nature.

**Hypothetical** (all Ring 1 reactions): The reaction is structurally plausible. Two compounds have the right structural relationship for this reaction type (e.g., they differ at exactly one stereocenter, so epimerization is possible). No database evidence has been checked.

**Inferred**: Not currently assigned by the pipeline. Reserved for reactions supported by indirect evidence (e.g., the D-form of this reaction is validated, so the L-form is inferred by mirror symmetry).

**Predicted**: The reaction was matched to a database entry with medium confidence. For example, a RHEA reaction was found where the substrate and product names fuzzy-match our compounds, but the match is not exact.

**Validated**: The reaction was matched to a database entry with high confidence. A RHEA reaction record directly confirms this substrate-product pair, or the match was manually reviewed and pinned.

Ring 2 enrichment is the mechanism that upgrades reactions from hypothetical to higher tiers. Without Ring 2, all reactions remain hypothetical.

## Data provenance

### Ring 1: generated from rules

All Ring 1 data is generated deterministically from stereochemistry rules. No external data is required.

**Compounds** are enumerated by iterating over all possible R/S configurations at each stereocenter. The number of stereocenters is determined by the carbon chain length and carbonyl position. This is exhaustive: every possible stereoisomer is generated.

**Names** come from `pipeline/data/name_mapping.json`, which maps stereo-configuration strings (like `"aldose-C6-RSSR"`) to known sugar names (like "D-Glucose"). Configurations not in the mapping get systematic names.

**Reactions** are generated by comparing all pairs (or in some cases all individual compounds) and checking structural criteria. For example, epimerization checks that two compounds of the same type and carbon count differ at exactly one stereocenter position. This is also exhaustive within each reaction type's domain.

**Formulas** are computed from carbon count and compound type:
- Aldoses/ketoses: C_n H_(2n) O_n
- Polyols: C_n H_(2n+2) O_n
- Each phosphate group adds: +1P, +3O, +1H (net effect of ester bond: H3PO4 minus H2O)

### Ring 2: external database enrichment

Ring 2 fetches data from four external databases and matches it to Ring 1 compounds and reactions.

**ChEBI** (Chemical Entities of Biological Interest): Provides compound identifiers, InChI strings, SMILES notation, and cross-references. SUGAR downloads the complete ChEBI dataset and matches compounds using a five-strategy engine:

1. Override pin: manually force a specific ChEBI ID (from `match_overrides.json`)
2. Override reject: manually block a match
3. Exact name: compound name exactly equals a ChEBI entry name
4. Alias match: one of the compound's aliases matches a ChEBI entry name
5. Formula unique: exactly one ChEBI entry has the same molecular formula
6. Fuzzy name: Levenshtein similarity ratio of at least 85%

Each match is assigned a confidence level (high, medium, or low). Match results are stored in `pipeline/cache/match_report.json` for review.

**KEGG** (Kyoto Encyclopedia of Genes and Genomes): Provides compound identifiers and pathway context. Fetched for each ChEBI-matched compound via the KEGG REST API.

**RHEA**: A curated database of biochemical reactions. SUGAR queries RHEA via SPARQL to find reactions involving matched compound pairs. When a RHEA reaction confirms a Ring 1 reaction, the evidence tier is upgraded to "validated" or "predicted" depending on match confidence.

**BRENDA**: A comprehensive enzyme database. SUGAR queries BRENDA's SOAP API to fetch EC numbers, enzyme names, and kinetic parameters for matched reactions. Requires authentication (email/password in `.env` file).

**D-to-L inference** (`import_/infer.py`): When a D-form reaction is validated by RHEA, the pipeline infers that the corresponding L-form reaction is also plausible. The L-form reaction gets the "inferred" evidence tier.

### Caching

All external API responses are cached locally in `pipeline/cache/`. This avoids repeated network requests and allows the pipeline to run offline after the first fetch. Cache can be selectively refreshed with `--refresh-chebi`, `--refresh-kegg`, `--refresh-rhea`, or `--refresh-brenda`, or fully refreshed with `--refresh`.

## Biochemistry rationale

### Stereochemistry model

SUGAR uses the CIP (Cahn-Ingold-Prelog) R/S system to describe stereocenters. Each monosaccharide is represented as a linear carbon chain with R or S at each chiral carbon.

**Aldoses** (carbonyl at C1): The first stereocenter is at C2. A C6 aldose has 4 stereocenters (C2, C3, C4, C5), producing 2^4 = 16 stereoisomers.

**Ketoses** (carbonyl at C2): The first stereocenter is at C3. A C6 ketose has 3 stereocenters (C3, C4, C5), producing 2^3 = 8 stereoisomers.

**D/L assignment**: Follows the Fischer convention. The configuration at the highest-numbered stereocenter determines D (R) or L (S). D-Glucose has stereocenters [R, S, S, R]; the last one is R, so it is D.

### Naming conventions

SUGAR uses a three-tier naming system:

1. **Known names**: Sugars with established common names (e.g., D-Glucose, D-Mannose, D-Fructose) are identified by their stereocenter pattern in `name_mapping.json`.
2. **Systematic names**: Unknown stereoisomers get names like "D-aldohexose-RSRS" indicating chirality, type, carbon count, and configuration.
3. **Phosphosugar names**: Derived from parent name with position suffix (e.g., "D-Glucose 6-phosphate"). Biologically important ones also get short aliases (e.g., "G6P").

### Reaction biochemistry

**Epimerization**: An epimerase enzyme inverts the configuration at a single stereocenter. The substrate and product differ at exactly one position in their stereocenter arrays. For example, D-Glucose [R,S,S,R] and D-Mannose [S,S,S,R] differ only at position 0 (C2). No cofactor required.

**Isomerization**: An isomerase converts an aldose to a ketose (or vice versa) by migrating the carbonyl from C1 to C2. This removes one stereocenter (the aldose C2 becomes the ketose carbonyl carbon). The remaining stereocenters must match. No cofactor required.

**Reduction**: An oxidoreductase (with NADH cofactor) reduces the carbonyl to a hydroxyl, producing a polyol. This can create a new stereocenter at the former carbonyl carbon. The reaction is modeled as irreversible because the reverse (oxidation) requires different conditions.

**Phosphorylation**: A kinase enzyme (with ATP cofactor) attaches a phosphate group to a specific hydroxyl position. Each position is treated independently. Biologically, C2 phosphorylation is excluded because that position is the carbonyl carbon.

**Dephosphorylation**: A phosphatase removes the phosphate group, regenerating the parent monosaccharide. No cofactor required.

**Mutase**: An intramolecular phosphotransferase moves a phosphate from one position to another on the same sugar. For example, glucose 1-phosphate to glucose 6-phosphate. Only applies between mono-phosphosugars that share the same parent and stereocenters.

**Phospho-epimerization and phospho-isomerization**: The same logic as standard epimerization and isomerization, but applied to phosphosugars. Both compounds must have the same phosphate positions. This prevents the generation of reactions that would require simultaneous stereocenter inversion and phosphate migration, which would not be a single enzymatic step.

### What is and is not included

**Included in Ring 1:**
- All possible C2-C7 aldose and ketose stereoisomers (exhaustive)
- All polyol reduction products with degeneracy handling
- Systematic C6 phosphosugars at positions 1, 3, 4, 6 (mono) and 1,6 / 3,6 (bis)
- Curated phosphosugars from other carbon lengths (important metabolic intermediates)
- All valid pairwise reactions within each reaction type

**Not included (some planned for future rings):**
- C2 phosphorylation (carbonyl position, not biologically relevant)
- C5 phosphorylation of hexoses (C5 hydroxyl is typically involved in ring closure)
- Sugar acids, lactones, amino sugars, nucleotide sugars, deoxy sugars (Ring 3)
- Disaccharides and oligosaccharides (Ring 4)
- Ring-form (pyranose/furanose) representations (SUGAR uses open-chain Fischer projection logic)
- Enzyme specificity constraints (SUGAR generates all structurally valid reactions, not just known enzyme-catalyzed ones)
- Thermodynamic favorability or kinetic rates (except what BRENDA provides in Ring 2)
- Multi-substrate reactions (all current reactions are single-substrate, single-product, ignoring cofactors as separate participants)

## Current coverage summary

| Category | Count | Details |
|----------|-------|---------|
| Aldoses | 63 | C2-C7, all stereoisomers |
| Ketoses | 31 | C3-C7, all stereoisomers |
| Polyols | 41 | Reduction products, degeneracy-resolved |
| Phosphosugars | 144 | 136 systematic + 8 curated |
| **Total compounds** | **279** | |
| Epimerizations | 478 | Bidirectional, same type and carbon count |
| Isomerizations | 124 | Bidirectional, aldose to ketose |
| Reductions | 94 | Irreversible, monosaccharide to polyol |
| Phosphorylations | 144 | One per phosphosugar |
| Dephosphorylations | 144 | One per phosphosugar |
| Mutases | 288 | Bidirectional, position migration |
| Phospho-epimerizations | 506 | Bidirectional, same positions |
| Phospho-isomerizations | 162 | Bidirectional, same positions |
| **Total reactions** | **1,940** | |
