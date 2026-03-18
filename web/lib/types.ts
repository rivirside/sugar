export type EvidenceTier = "validated" | "predicted" | "inferred" | "hypothetical";
export type CompoundType =
  | "aldose"
  | "ketose"
  | "polyol"
  | "phosphate"
  | "acid"
  | "lactone"
  | "amino_sugar"
  | "nucleotide_sugar"
  | "deoxy_sugar"
  | "disaccharide";
export type Chirality = "D" | "L" | "achiral";
export type ReactionType =
  | "epimerization"
  | "isomerization"
  | "reduction"
  | "oxidation"
  | "phosphorylation"
  | "dephosphorylation"
  | "mutase"
  | "lactonization"
  | "aldol"
  | "hydrolysis"
  | "condensation"
  | "transamination"
  | "nucleotidyltransfer"
  | "multi_epimerization";

export interface Compound {
  id: string;
  name: string;
  aliases: string[];
  type: CompoundType;
  carbons: number;
  chirality: Chirality;
  formula: string;
  stereocenters: string[];
  modifications: Record<string, unknown> | null;
  parent_monosaccharide: string | null;
  commercial: boolean;
  cost_usd_per_kg: number | null;
  metadata: Record<string, unknown>;
}

export interface Reaction {
  id: string;
  substrates: string[];
  products: string[];
  reaction_type: ReactionType;
  evidence_tier: EvidenceTier;
  evidence_criteria: Record<string, unknown> | unknown[];
  yield: number | null;
  cofactor_burden: number;
  cost_score: number;
  // Optional fields present in richer datasets
  ec_number?: string | null;
  enzyme_id?: string | null;
  enzyme_name?: string | null;
  cofactors?: string[];
  delta_g?: number | null;
  conditions?: Record<string, unknown> | null;
  pmid?: string[];
  rhea_id?: string | null;
  organism?: string[];
  km_mm?: number | null;
  kcat_sec?: number | null;
  metadata?: Record<string, unknown>;
}

export interface PipelineMetadata {
  pipeline_version: string;
  generated_at: string;
  counts: {
    total_compounds: number;
    total_reactions: number;
    monosaccharides?: number;
    polyols?: number;
    epimerizations?: number;
    isomerizations?: number;
    reductions?: number;
    [key: string]: number | undefined;
  };
  completeness_warnings?: unknown[];
  duplicate_warnings?: number;
  output_files?: Record<string, string>;
}

export interface PathwayStep {
  reaction: Reaction;
  from_compound: Compound;
  to_compound: Compound;
}

export interface Pathway {
  steps: PathwayStep[];
  total_cost: number;
  total_yield: number;
  step_count: number;
}
