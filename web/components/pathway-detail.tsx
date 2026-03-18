"use client";

import Link from "next/link";
import { EvidenceBadge } from "@/components/evidence-badge";
import { compoundMap } from "@/lib/data";
import { formatYield, cumulativeYield } from "@/lib/utils";
import type { PathResult } from "@/lib/pathfinding";
import type { Reaction } from "@/lib/types";
import { ArrowRight } from "lucide-react";

interface PathwayDetailProps {
  pathway: PathResult;
  reactionMap: Map<string, Reaction>;
}

export function PathwayDetail({ pathway, reactionMap }: PathwayDetailProps) {
  const reactions = pathway.reactionIds
    .map((id) => reactionMap.get(id))
    .filter(Boolean) as Reaction[];
  const yields = reactions.map((r) => r.yield);
  const cumYield = cumulativeYield(yields);

  return (
    <div className="flex flex-col gap-4">
      <h3 className="text-sm font-medium text-zinc-300">
        Pathway Detail ({pathway.nodes.length - 1} steps)
      </h3>

      <div className="flex flex-col gap-3">
        {reactions.map((reaction, i) => {
          const fromId = pathway.nodes[i];
          const toId = pathway.nodes[i + 1];
          const from = compoundMap.get(fromId);
          const to = compoundMap.get(toId);

          return (
            <div
              key={reaction.id}
              className="rounded-lg border border-zinc-800 bg-zinc-900/60 p-3"
            >
              <div className="flex items-center gap-2 text-xs text-zinc-500">
                <span className="font-medium">Step {i + 1}</span>
                <span className="capitalize">{reaction.reaction_type.replace("_", " ")}</span>
                <EvidenceBadge tier={reaction.evidence_tier} />
              </div>

              <div className="mt-2 flex items-center gap-2">
                <Link
                  href={`/compound/${fromId}`}
                  className="text-sm font-medium text-zinc-200 hover:text-white"
                >
                  {from?.name ?? fromId}
                </Link>
                <ArrowRight className="h-3.5 w-3.5 text-zinc-600" />
                <Link
                  href={`/compound/${toId}`}
                  className="text-sm font-medium text-zinc-200 hover:text-white"
                >
                  {to?.name ?? toId}
                </Link>
              </div>

              <div className="mt-2 flex gap-4 text-xs text-zinc-500">
                <span>Yield: {formatYield(reaction.yield)}</span>
                <span>Cost: {reaction.cost_score.toFixed(2)}</span>
                {reaction.cofactors && reaction.cofactors.length > 0 && (
                  <span>Cofactors: {reaction.cofactors.join(", ")}</span>
                )}
              </div>
            </div>
          );
        })}
      </div>

      {/* Summary */}
      <div className="rounded-lg border border-zinc-800 bg-zinc-800/40 p-3">
        <h4 className="text-xs font-medium uppercase tracking-wider text-zinc-500">
          Summary
        </h4>
        <div className="mt-2 grid grid-cols-3 gap-4 text-sm">
          <div>
            <span className="text-zinc-500">Steps</span>
            <p className="font-medium text-zinc-200">
              {pathway.nodes.length - 1}
            </p>
          </div>
          <div>
            <span className="text-zinc-500">Est. Yield</span>
            <p className="font-medium text-zinc-200">{formatYield(cumYield)}</p>
          </div>
          <div>
            <span className="text-zinc-500">Total Cost</span>
            <p className="font-medium text-zinc-200">
              {pathway.totalCost.toFixed(2)}
            </p>
          </div>
        </div>
      </div>
    </div>
  );
}
