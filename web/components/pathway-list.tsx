"use client";

import { EvidenceBadge } from "@/components/evidence-badge";
import { compoundMap } from "@/lib/data";
import { formatYield, cumulativeYield } from "@/lib/utils";
import type { PathResult } from "@/lib/pathfinding";
import type { Reaction } from "@/lib/types";

interface PathwayListProps {
  pathways: PathResult[];
  reactionMap: Map<string, Reaction>;
  selectedIndex: number;
  onSelect: (index: number) => void;
}

export function PathwayList({
  pathways,
  reactionMap,
  selectedIndex,
  onSelect,
}: PathwayListProps) {
  if (pathways.length === 0) {
    return (
      <div className="flex h-48 items-center justify-center text-sm text-zinc-500">
        No pathways found
      </div>
    );
  }

  return (
    <div className="flex flex-col gap-2">
      {pathways.map((pathway, index) => {
        const reactions = pathway.reactionIds
          .map((id) => reactionMap.get(id))
          .filter(Boolean) as Reaction[];
        const yields = reactions.map((r) => r.yield);
        const cumYield = cumulativeYield(yields);
        const worstTier = reactions.reduce<string>((worst, r) => {
          const order = ["validated", "predicted", "inferred", "hypothetical"];
          return order.indexOf(r.evidence_tier) > order.indexOf(worst)
            ? r.evidence_tier
            : worst;
        }, "validated");

        // Build condensed chain
        const chain = pathway.nodes
          .map((id) => compoundMap.get(id)?.name ?? id)
          .join(" \u2192 ");

        return (
          <button
            key={index}
            type="button"
            onClick={() => onSelect(index)}
            className={`w-full rounded-lg border p-3 text-left transition-colors ${
              index === selectedIndex
                ? "border-zinc-600 bg-zinc-800/80"
                : "border-zinc-800 bg-zinc-900/60 hover:border-zinc-700"
            }`}
          >
            <div className="flex items-center gap-2">
              <span className="flex h-5 w-5 shrink-0 items-center justify-center rounded bg-zinc-800 text-xs font-medium text-zinc-400">
                {index + 1}
              </span>
              <span className="text-xs text-zinc-400">
                {pathway.nodes.length - 1} step
                {pathway.nodes.length - 1 !== 1 ? "s" : ""}
              </span>
              <EvidenceBadge
                tier={worstTier as "validated" | "predicted" | "inferred" | "hypothetical"}
              />
              <span className="ml-auto text-xs text-zinc-500">
                yield: {formatYield(cumYield)}
              </span>
            </div>
            <p className="mt-2 truncate text-xs text-zinc-300">{chain}</p>
          </button>
        );
      })}
    </div>
  );
}
