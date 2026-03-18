import compoundsJson from "@/data/compounds.json";
import reactionsJson from "@/data/reactions.json";
import metadataJson from "@/data/pipeline_metadata.json";
import type { Compound, Reaction, PipelineMetadata } from "./types";

export const compounds: Compound[] = compoundsJson as Compound[];
export const reactions: Reaction[] = reactionsJson as Reaction[];
export const metadata: PipelineMetadata = metadataJson as PipelineMetadata;

export const compoundMap = new Map<string, Compound>(
  compounds.map((c) => [c.id, c])
);
export const reactionMap = new Map<string, Reaction>(
  reactions.map((r) => [r.id, r])
);

export function getReactionsForCompound(compoundId: string): Reaction[] {
  return reactions.filter(
    (r) =>
      r.substrates.includes(compoundId) || r.products.includes(compoundId)
  );
}
