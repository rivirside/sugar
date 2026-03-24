import type { Reaction, ScoringMode } from "./types";

export interface Edge {
  reactionId: string;
  target: string;
  weight: number;
}

export type AdjacencyList = Map<string, Edge[]>;

function getWeight(r: Reaction, mode: ScoringMode, alpha: number): number {
  if (mode === "cost") return r.cost_score;
  if (mode === "engineerability") return r.engineerability_score ?? r.cost_score;
  // combined
  const eng = r.engineerability_score ?? r.cost_score;
  return alpha * r.cost_score + (1 - alpha) * eng;
}

export function buildGraph(
  reactions: Reaction[],
  scoringMode: ScoringMode = "cost",
  alpha: number = 0.5,
): AdjacencyList {
  const adj: AdjacencyList = new Map();
  for (const r of reactions) {
    if (r.substrates.length !== 1 || r.products.length !== 1) continue;
    const source = r.substrates[0];
    const target = r.products[0];
    const weight = getWeight(r, scoringMode, alpha);
    if (!adj.has(source)) adj.set(source, []);
    adj.get(source)!.push({ reactionId: r.id, target, weight });
  }
  return adj;
}
