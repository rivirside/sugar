# Enzyme Gap Analysis Design Spec (Ring 4)

## Overview

Add an enzyme gap analysis stage to the SUGAR pipeline that evaluates every reaction for enzyme availability and engineering feasibility. This transforms SUGAR from a pathway encyclopedia into a protein engineering decision tool.

The primary use case: finding the most buildable multi-step route between two compounds (e.g., D-Glucose to L-Glucose) by ranking routes not just by biochemical cost, but by how many steps have known enzymes and how engineerable the gaps are.

## Goals

1. Classify every reaction's enzyme coverage level (direct match, cross-substrate candidate, EC family only, none)
2. For reactions without direct enzyme matches, find structurally similar enzymes that could be engineering starting points
3. Score each reaction's engineerability on a 0.0-1.0 scale
4. Enable pathfinding to rank routes by engineering feasibility
5. Cache a lightweight enzyme index for fast route comparison; defer full dossier data to on-demand fetching

## Non-Goals (for this phase)

- Frontend workbench UI (will be built later on top of this data layer)
- Live UniProt/PDB sequence/structure fetching (deferred to on-demand layer)
- Directed evolution strategy recommendations (out of scope)
- Enzyme kinetics optimization (Ring 2 BRENDA data covers this separately)

## Pipeline Placement

Ring 4 runs after Ring 2 enrichment and Ring 3 phosphosugars. It reads all compounds and reactions (including any Ring 2 enzyme annotations) and adds new fields to each reaction. It does not create or remove compounds or reactions.

```
Ring 1: Enumerate + React (compounds, reactions)
Ring 2: Import + Enrich (external database annotations)
Ring 3: Phosphosugars (derivative compounds + reactions)
Ring 4: Enzyme Gap Analysis (coverage classification, cross-substrate matching, engineerability scoring)  <-- NEW
Validate + Export
```

Ring 4 depends on Ring 2 data (enzyme annotations from RHEA/BRENDA). If Ring 2 was skipped, all reactions will be classified as "none" or "family_only" depending on EC family data availability.

## Cross-Substrate Matching

### Problem

When a reaction has no known enzyme (evidence_tier = "hypothetical"), we need to find the closest existing enzyme that could be engineered to catalyze it. SUGAR already knows the structural relationship between every pair of compounds, so it can quantify how similar a known enzyme's substrate is to the target substrate.

### Matching Layers

For a gap reaction (e.g., "epimerize position C3 on D-Allose-6-phosphate"), candidates are searched in three layers ordered by relevance:

**Layer 1: Same reaction type + same position.**
Find all validated/predicted reactions that invert the same stereocenter at the same position. Example: a known C3 epimerase on D-Glucose is a Layer 1 candidate for a C3 epimerization gap on D-Allose. These are the best candidates because the enzyme already performs the right chemistry at the right position.

**Layer 2: Same reaction type + different position.**
Find all validated/predicted reactions of the same type but at a different position. A known C2 epimerase is still an epimerase, but redirecting it to C3 requires more engineering. Worse than Layer 1, but still relevant.

**Layer 3: Same EC subclass.**
Find all enzymes in the same EC subclass (e.g., all 5.1.3.x epimerases) from BRENDA/UniProt, even if SUGAR does not model the exact reaction. This catches enzymes that work on substrates outside SUGAR's current scope.

### Position Extraction

Position extraction depends on reaction type:

**Epimerization**: compare substrate and product stereocenters arrays. Exactly one position differs. That index is the position.

**Mutase**: stereocenters are identical between substrate and product (phosphate moves, chirality doesn't change). Position is extracted by comparing `modifications[].position` on substrate vs product. The "from" and "to" positions together define the mutase position pair.

**Phosphorylation / Dephosphorylation**: position is the phosphate site from `modifications[].position` on the phosphorylated compound. A C1 kinase and a C6 kinase operate at different positions and should be distinguished.

**Isomerization**: no position concept (aldose-ketose conversion involves the carbonyl, not a specific stereocenter). Layers 1 and 2 collapse into a single "same reaction type" layer.

**Reduction**: no position concept (carbonyl reduction). Layers 1 and 2 collapse into a single "same reaction type" layer.

### Reaction Type Mapping

The pipeline stores phospho-epimerizations as `reaction_type: "epimerization"` and phospho-isomerizations as `reaction_type: "isomerization"` (not distinct types). This means cross-substrate matching for an epimerization gap will include both base-sugar and phosphosugar epimerases in the candidate pool. The substrate similarity score's modification_distance dimension naturally penalizes candidates whose substrates have different phosphate patterns, so phosphosugar candidates will rank lower for base-sugar gaps and vice versa. No additional filtering is needed.

Note: the pipeline generates 6 distinct `reaction_type` values: "epimerization", "isomerization", "reduction", "phosphorylation", "dephosphorylation", "mutase". The 8 reaction generators map to these 6 types (phospho-epimerization -> "epimerization", phospho-isomerization -> "isomerization").

### Directionality

Epimerization is inherently reversible (enzyme that epimerizes A to B also epimerizes B to A). Mutase is also reversible (same enzyme, different phosphate position equilibrium). For these symmetric reaction types, directionality is ignored in matching.

For directional reaction types (phosphorylation vs dephosphorylation), directionality is enforced. A kinase candidate is not relevant for a phosphatase gap.

Note on reduction: the pipeline generates reduction reactions but not oxidation reactions. The reverse direction does not exist in the dataset. Cross-substrate matching for a reduction gap will find other reduction reactions as candidates. If a future use case needs oxidation matching, oxidation reactions would need to be added to the pipeline first.

## Multi-Dimensional Substrate Similarity

### Problem

A single-number similarity score that only counts stereocenter differences would be misleading. An enzyme that works on D-Glucose is a better engineering starting point for D-Allose (same carbon count, same type, just stereocenters differ) than an enzyme that works on D-Glucose-6-phosphate (bulky charged phosphate group the active site was shaped around).

### Dimensions

**Dimension 1: Stereocenter distance.**
Count of stereocenters that differ between the known enzyme's substrate and the target substrate. This is the easiest kind of difference to engineer around (reshaping a binding pocket). Range: 0 to max(len(stereocenters_a), len(stereocenters_b)).

**Dimension 2: Modification distance.**
Differences in phosphate groups (or future modification types). Scored as:
- Same modifications (identical positions and types): 0.0
- Same type, different position (e.g., C1-phosphate vs C6-phosphate): 0.4
- Different number of modifications (has phosphate vs doesn't): 0.7
- Different modification types (future-proofing): 0.9
Range: 0.0 to 1.0.

**Dimension 3: Carbon count distance.**
Absolute difference in chain length divided by max chain length for normalization. An enzyme shaped for C6 may tolerate C5 but C3 is a stretch. Range: 0.0 to 1.0.

**Dimension 4: Type distance.**
Aldose vs ketose vs polyol. Same type = 0.0. Aldose-ketose = 0.5 (carbonyl position changes, but carbon skeleton similar). Either to polyol = 0.7 (reduced carbonyl, different chemistry). Range: 0.0 to 1.0.

### Composite Score

```
# Normalization:
#   stereo_norm = stereo_distance / max_stereo if max_stereo > 0 else 0.0
#   carbon_norm = abs(carbons_a - carbons_b) / max(carbons_a, carbons_b)
# Edge case: when both compounds have 0 stereocenters, stereo_norm = 0.0 (identical)

substrate_similarity = 1.0 - (
    W_stereo * stereo_norm +
    W_mod    * modification_distance +
    W_carbon * carbon_norm +
    W_type   * type_distance
)
```

Weights:
- W_stereo = 0.35 (stereocenters are the most common and most engineerable difference)
- W_mod = 0.25 (modifications significantly change binding pocket requirements)
- W_carbon = 0.25 (chain length affects active site geometry)
- W_type = 0.15 (type differences are important but less common in cross-substrate searches)

Score range: 0.0 (maximally different) to 1.0 (identical substrate).

The per-dimension breakdown is stored alongside the composite score so users can judge which dimensions matter most for their engineering approach.

## Engineerability Score

### Components

**Component 1: Coverage level (weight: 0.4)**
Base penalty by coverage level, modulated by candidate count:
- Direct enzyme match: 0.0
- Cross-substrate Layer 1 candidate exists: base 0.3
- Cross-substrate Layer 2 candidate exists: base 0.6
- Layer 3 only (EC family, no substrate data): base 0.8
- No candidates at all: 1.0

When candidates exist (Layers 1-3), the base penalty is reduced by up to 30% based on candidate count: `adjusted = base * (1.0 - min(num_candidates / 5, 1.0) * 0.3)`. Having 5+ candidates at Layer 1 reduces the penalty from 0.3 to 0.21 (more starting points = easier engineering).

**Component 2: Best candidate substrate similarity (weight: 0.3)**
The multi-dimensional similarity score of the top-ranked candidate, inverted (1.0 - similarity) so that higher = harder. Set to 1.0 when no candidates exist. Range: 0.0 to 1.0.

**Component 3: EC family richness (weight: 0.15)**
Number of distinct enzymes in the relevant EC subclass, normalized. Large family = more starting points for engineering. Computed as: 1.0 - min(family_size / 50, 1.0). A family of 50+ enzymes scores 0.0 (easy). A family of 1 scores 0.98. No family scores 1.0. Range: 0.0 to 1.0.

**Component 4: Structural data availability (weight: 0.15)**
Does the best candidate have a crystal structure in PDB? Yes = 0.0. No = 1.0. Range: 0.0 or 1.0 (binary in phase 1). AlphaFold model availability may be added in a future phase as an intermediate tier (0.3), but phase 1 only checks PDB since that data is available from the enzyme index.

### Formula

```
engineerability_score = (
    0.4 * coverage_level +
    0.3 * (1.0 - best_similarity) +
    0.15 * family_richness_penalty +
    0.15 * structural_data_penalty
)
```

Range: 0.0 (trivially engineerable, enzyme exists) to 1.0 (no leads at all).

## Integration with Pathfinding

### Separation from cost_score

The engineerability_score is kept separate from the existing cost_score. This avoids breaking existing behavior and enables distinct ranking modes:

- **Biochemical cost** (current behavior): rank by yield, cofactor burden, evidence tier
- **Engineerability** (new): rank by how buildable the route is given known enzymes
- **Combined** (new): weighted blend of both scores

The pathfinder gains a `scoring_mode` parameter: "cost" (default, backward-compatible), "engineerability", or "combined".

### Combined scoring

```
combined_score = alpha * cost_score + (1 - alpha) * engineerability_score
```

Default alpha = 0.5. Configurable by the user. When alpha = 1.0, behaves identically to current cost-only mode.

## Data Model Changes

### Reaction fields (new)

Each reaction dict gains:

```python
{
    # Existing fields unchanged...

    # New Ring 4 fields:
    "enzyme_coverage": str,           # "direct" | "cross_substrate" | "family_only" | "none"
                                      # Classification rule: "direct" if ec_number present from Ring 2,
                                      # "cross_substrate" if best candidate is Layer 1 or 2,
                                      # "family_only" if best candidate is Layer 3 only,
                                      # "none" if no candidates found

    "cross_substrate_candidates": [   # Sorted by similarity descending, top 5
        {
            "ec_number": str,
            "enzyme_name": str,
            "organism": str | None,
            "uniprot_id": str | None,
            "pdb_ids": list[str],
            "source_reaction_id": str,
            "known_substrate_id": str,
            "matching_layer": int,     # 1, 2, or 3
            "similarity": {
                "overall": float,       # 0.0-1.0
                "stereocenter_distance": int,
                "modification_distance": float,
                "carbon_count_distance": int,
                "type_distance": float
            }
        }
    ],

    "ec_family_size": int | None,

    "engineerability_score": float,    # 0.0-1.0

    "engineerability_components": {
        "coverage_level": float,
        "best_similarity": float,
        "family_richness": float,
        "structural_data": float
    }
}
```

### New output file: enzyme_index.json

Lightweight enzyme family index keyed by EC number. Built in two tiers:

**Tier 1 (from Ring 2 data, always available):**
```python
{
    "5.1.3.2": {
        "name": "UDP-glucose 4-epimerase",      # From Ring 2 enzyme_name
        "organisms": ["E. coli", "H. sapiens"],  # From Ring 2 reaction organisms
        "known_substrates": ["D-GLC", "D-GAL"],  # From Ring 2 reaction substrates
        "reaction_count": 3                       # Number of SUGAR reactions with this EC
    }
}
```

**Tier 2 (from BRENDA/UniProt API, available when credentials configured):**
Adds to the Tier 1 entry:
```python
{
    "family_size": 47,                           # Distinct enzymes in EC subclass
    "pdb_count": 12,                             # PDB structures available
    "uniprot_ids": ["P09147", ...]               # Representative entries only
}
```

When BRENDA credentials are unavailable, Tier 2 fields are set to None. The engineerability score's Component 3 (family richness) defaults to 0.8 and Component 4 (structural data) defaults to 1.0 in this case. This produces conservative (high) engineerability scores, which is correct: less data means less certainty about engineering feasibility.

This index supports the summary view (green/yellow/red per step). The full dossier (sequences, structures, variant lists) is fetched on demand by the frontend in a future phase.

### pipeline_metadata.json additions

```python
{
    # Existing fields...
    "gap_analysis": {
        "reactions_analyzed": int,
        "coverage_direct": int,
        "coverage_cross_substrate": int,
        "coverage_family_only": int,
        "coverage_none": int,
        "avg_engineerability_score": float,
        "ec_families_indexed": int
    }
}
```

## Module Structure

```
pipeline/
  analyze/
    __init__.py
    gap_analysis.py      # Ring 4 orchestrator: classify, match, score
    cross_substrate.py   # Layer 1-3 matching logic
    similarity.py        # Multi-dimensional substrate similarity scoring
    engineerability.py   # Composite score computation
    enzyme_index.py      # EC family metadata fetching and caching
```

### gap_analysis.py

Main entry point: `run_gap_analysis(compounds, reactions, enzyme_index=None) -> tuple[list[dict], dict]`

Returns updated reactions list and gap_analysis metadata dict.

Steps:
1. Build compound lookup by ID
2. Build reaction-type index (group reactions by reaction_type)
3. For each reaction, classify enzyme_coverage based on existing Ring 2 fields
4. For reactions without direct coverage, call cross_substrate.find_candidates()
5. Compute engineerability_score via engineerability.compute_score()
6. Return enriched reactions and summary metadata

### cross_substrate.py

Main entry point: `find_candidates(gap_reaction, all_reactions, compound_map, max_candidates=5) -> list[dict]`

Logic:
1. Determine the gap reaction's type, position (if applicable), and substrate compound
2. Layer 1: filter all_reactions for same type + same position + has enzyme data. Score each by substrate similarity.
3. Layer 2: filter for same type + different position + has enzyme data. Score each.
4. Layer 3: look up EC subclass in enzyme_index for broader family matches.
5. Merge, deduplicate by EC number (keeping the entry with the best rank: lowest layer, then highest similarity), sort by (layer ascending, similarity descending), take top max_candidates.

"Has enzyme data" is defined as: `ec_number is not None`. This is the broadest useful filter. A reaction with an EC number but no enzyme_name still represents a known enzymatic activity.

### similarity.py

Main entry point: `compute_similarity(compound_a, compound_b) -> dict`

Returns:
```python
{
    "overall": float,
    "stereocenter_distance": int,
    "modification_distance": float,
    "carbon_count_distance": int,
    "type_distance": float
}
```

Pure function, no side effects. Operates on compound dicts.

### engineerability.py

Main entry point: `compute_score(coverage_level: str, best_similarity: float, ec_family_size: int | None, has_pdb: bool, num_candidates: int = 0) -> tuple[float, dict]`

When `ec_family_size` is None (no BRENDA data available), Component 3 defaults to 0.8 (assumes small family, penalized but not maximally). The `num_candidates` parameter is used by Component 1 to modulate the base penalty by candidate count.

Returns (engineerability_score, components_dict).

Pure function, no side effects.

### enzyme_index.py

Main entry point: `build_enzyme_index(reactions) -> dict`

Builds the EC-number-keyed index from reactions that have Ring 2 enzyme annotations. If BRENDA/UniProt API data is available in the cache, incorporates family size and PDB counts.

Secondary entry point: `fetch_ec_family_data(ec_number, cache_dir) -> dict`

Fetches family metadata from BRENDA/UniProt for a given EC number. Caches responses to avoid repeated API calls.

## Pipeline Orchestrator Changes

`run_pipeline.py` gains new steps after Ring 3:

```
[G1] Build enzyme index from Ring 2 data
[G2] Run gap analysis (classify, match, score)
[G3] Export enzyme_index.json
```

Step numbering follows the Ring 2 convention with a prefix (`[G1]`, `[G2]`, `[G3]` for Gap analysis) to avoid renumbering existing steps.

## Pathfinder Changes

### pipeline/reactions/score.py

Add engineerability-aware scoring:

```python
def compute_combined_score(cost_score, engineerability_score, alpha=0.5):
    return alpha * cost_score + (1 - alpha) * engineerability_score
```

### web/lib/graph.ts

The graph builder gains support for scoring_mode:
- "cost": use cost_score as edge weight (current behavior, default)
- "engineerability": use engineerability_score as edge weight
- "combined": use combined_score with configurable alpha

The pathfinding functions (Dijkstra, Yen's) remain unchanged. They operate on edge weights, and the graph builder selects which score to use as the weight.

### Multi-substrate reactions

The current `buildGraph()` in graph.ts skips reactions where `substrates.length !== 1 || products.length !== 1`. These reactions cannot appear in pathfinding results. Ring 4 still computes engineerability scores for multi-substrate reactions (the data is useful for the workbench view), but they are excluded from route ranking.

## TypeScript Type Changes

The `Reaction` interface in `web/lib/types.ts` gains optional Ring 4 fields:

```typescript
// Ring 4: Enzyme Gap Analysis (optional, present when Ring 4 has run)
enzyme_coverage?: "direct" | "cross_substrate" | "family_only" | "none";
cross_substrate_candidates?: Array<{
  ec_number: string;
  enzyme_name: string;
  organism: string | null;
  uniprot_id: string | null;
  pdb_ids: string[];
  source_reaction_id: string;
  known_substrate_id: string;
  matching_layer: 1 | 2 | 3;
  similarity: {
    overall: number;
    stereocenter_distance: number;
    modification_distance: number;
    carbon_count_distance: number;
    type_distance: number;
  };
}>;
ec_family_size?: number | null;
engineerability_score?: number;
engineerability_components?: {
  coverage_level: number;
  best_similarity: number;
  family_richness: number;
  structural_data: number;
};
```

The `ScoringMode` type is added: `type ScoringMode = "cost" | "engineerability" | "combined";`

The `buildGraph` function signature gains an optional `scoringMode` parameter.

## Testing Strategy

### Unit tests: similarity.py

- Identical compounds: all dimensions 0, overall 1.0
- One stereocenter different, everything else same: high similarity
- Different carbon count: carbon_count_distance > 0, lower overall
- Same sugar with vs without phosphate: modification_distance > 0
- Aldose vs ketose: type_distance = 0.5
- Maximally different compounds (different type, carbon count, stereocenters, modifications): low overall
- Both compounds have zero stereocenters (e.g., Glycolaldehyde C2): stereo_norm = 0.0, no division by zero

### Unit tests: cross_substrate.py

- Gap with direct enzyme match: classified as "direct," no cross-substrate search
- Gap with Layer 1 match (same type, same position, different substrate): found and ranked
- Gap with Layer 2 match only (same type, different position): found at Layer 2
- Gap with Layer 3 match only (EC family, no substrate data): found at Layer 3
- Gap with no match: empty candidates, coverage "none"
- Multiple candidates: sorted by (layer asc, similarity desc), capped at max_candidates

### Unit tests: engineerability.py

- Direct enzyme: score 0.0
- Cross-substrate Layer 1 with high similarity and large family and PDB: low score
- No candidates at all: score 1.0
- Boundary: all components at maximum penalty = 1.0
- Boundary: all components at zero = 0.0

### Integration test

- Run full pipeline with Ring 4 enabled on existing dataset
- Every reaction has enzyme_coverage field
- Reactions with Ring 2 enzyme data classified as "direct"
- Reactions without classified and scored
- enzyme_index.json is valid JSON with expected structure
- pipeline_metadata.json contains gap_analysis section
- D-GLC to L-GLC pathfinding returns routes ranked by engineerability
- Engineerability ranking differs from cost ranking (proving the new score adds information)

### Unit tests: polyol and phosphosugar edge cases

- Polyol-to-polyol path: all reduction steps score near 1.0 engineerability (no known reductases in dataset). Combined scoring falls back to cost ranking.
- Phosphosugar epimerization gap: cross-substrate matching finds both base-sugar and phosphosugar epimerases. Candidates with matching modification pattern rank higher.
- Mutase gap: position extracted from modifications, not stereocenters. Candidates matched by phosphate position pair.

### Offline testing

All tests run offline with mocked external API responses. No live BRENDA/UniProt calls in the test suite. Mock data includes realistic EC family sizes and PDB counts for test coverage of all scoring branches.

## Open Questions (Resolved)

**Q: Should Ring 4 run if Ring 2 was skipped?**
A: Yes. All reactions will be classified as "none" or "family_only." The engineerability scores will all be high (near 1.0), which is correct: without enzyme data, everything is hard to engineer. The pathfinder's "engineerability" mode will produce flat rankings, naturally falling back to cost-based ranking when combined.

**Q: How many cross-substrate candidates to store per reaction?**
A: Top 5, sorted by (layer ascending, similarity descending). This balances data size with usefulness. The full candidate list can be recomputed on demand if needed.

**Q: Should the similarity weights be configurable?**
A: Not in phase 1. Hardcode the weights (W_stereo=0.35, W_mod=0.25, W_carbon=0.25, W_type=0.15). If users need to tune them, add a config parameter later.
