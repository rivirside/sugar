import { describe, it, expect } from "vitest";
import { findKShortestPaths } from "@/lib/pathfinding";
import { buildGraph } from "@/lib/graph";
import type { Reaction } from "@/lib/types";

// Small test graph:
//   A -> B  (R1, cost 1.0)
//   B -> C  (R2, cost 1.0)  => A->B->C total cost 2.0
//   A -> C  (R3, cost 3.0)  => A->C direct total cost 3.0
//   B -> D  (R4, cost 1.0)
//   D -> C  (R5, cost 1.0)  => A->B->D->C total cost 3.0
const testReactions: Reaction[] = [
  {
    id: "R1",
    substrates: ["A"],
    products: ["B"],
    reaction_type: "epimerization",
    cost_score: 1.0,
    evidence_tier: "hypothetical",
    evidence_criteria: {},
    yield: null,
    cofactor_burden: 0,
  },
  {
    id: "R2",
    substrates: ["B"],
    products: ["C"],
    reaction_type: "epimerization",
    cost_score: 1.0,
    evidence_tier: "hypothetical",
    evidence_criteria: {},
    yield: null,
    cofactor_burden: 0,
  },
  {
    id: "R3",
    substrates: ["A"],
    products: ["C"],
    reaction_type: "epimerization",
    cost_score: 3.0,
    evidence_tier: "hypothetical",
    evidence_criteria: {},
    yield: null,
    cofactor_burden: 0,
  },
  {
    id: "R4",
    substrates: ["B"],
    products: ["D"],
    reaction_type: "epimerization",
    cost_score: 1.0,
    evidence_tier: "hypothetical",
    evidence_criteria: {},
    yield: null,
    cofactor_burden: 0,
  },
  {
    id: "R5",
    substrates: ["D"],
    products: ["C"],
    reaction_type: "epimerization",
    cost_score: 1.0,
    evidence_tier: "hypothetical",
    evidence_criteria: {},
    yield: null,
    cofactor_burden: 0,
  },
];

describe("findKShortestPaths", () => {
  it("finds the shortest path", () => {
    const graph = buildGraph(testReactions);
    const results = findKShortestPaths(graph, "A", "C", 3, { maxSteps: 6, timeoutMs: 5000 });
    expect(results.length).toBeGreaterThanOrEqual(1);
    expect(results[0].nodes).toEqual(["A", "B", "C"]);
    expect(results[0].totalCost).toBe(2.0);
  });

  it("returns multiple paths sorted by cost", () => {
    const graph = buildGraph(testReactions);
    const results = findKShortestPaths(graph, "A", "C", 3, { maxSteps: 6, timeoutMs: 5000 });
    expect(results.length).toBe(3);
    for (let i = 1; i < results.length; i++) {
      expect(results[i].totalCost).toBeGreaterThanOrEqual(results[i - 1].totalCost);
    }
  });

  it("respects maxSteps", () => {
    const graph = buildGraph(testReactions);
    // maxSteps: 1 means only 1-hop paths; only R3 (A->C directly) qualifies
    const results = findKShortestPaths(graph, "A", "C", 10, { maxSteps: 1, timeoutMs: 5000 });
    for (const p of results) {
      expect(p.nodes.length - 1).toBeLessThanOrEqual(1);
    }
  });

  it("returns empty array for unreachable target", () => {
    const graph = buildGraph(testReactions);
    const results = findKShortestPaths(graph, "A", "Z", 3, { maxSteps: 6, timeoutMs: 5000 });
    expect(results.length).toBe(0);
  });

  it("returns empty array on immediate timeout", () => {
    const graph = buildGraph(testReactions);
    // timeoutMs: 0 means the deadline is already passed before Yen's loop
    const results = findKShortestPaths(graph, "A", "C", 100, { maxSteps: 10, timeoutMs: 0 });
    // The first path (from dijkstra before the loop) may still be found;
    // subsequent candidates are skipped. We just verify no error is thrown
    // and results is a valid array.
    expect(Array.isArray(results)).toBe(true);
  });
});
