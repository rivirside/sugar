"use client";

import { useState, useEffect, useMemo, useCallback, Suspense } from "react";
import { useRouter, useSearchParams } from "next/navigation";
import { CompoundSearch } from "@/components/compound-search";
import { PathwayList } from "@/components/pathway-list";
import { PathwayDetail } from "@/components/pathway-detail";
import { reactions, compoundMap, reactionMap } from "@/lib/data";
import { buildGraph } from "@/lib/graph";
import { findKShortestPaths, type PathResult } from "@/lib/pathfinding";
import type { Compound, EvidenceTier } from "@/lib/types";
import { ArrowRight, Loader2 } from "lucide-react";

const EVIDENCE_TIERS: EvidenceTier[] = [
  "validated",
  "predicted",
  "inferred",
  "hypothetical",
];

function PathwayFinderContent() {
  const router = useRouter();
  const searchParams = useSearchParams();

  const sourceId = searchParams.get("source") ?? "";
  const targetId = searchParams.get("target") ?? "";

  const [source, setSource] = useState<Compound | null>(
    sourceId ? compoundMap.get(sourceId) ?? null : null
  );
  const [target, setTarget] = useState<Compound | null>(
    targetId ? compoundMap.get(targetId) ?? null : null
  );

  const [maxSteps, setMaxSteps] = useState(6);
  const [enabledTiers, setEnabledTiers] = useState<Set<EvidenceTier>>(
    new Set(EVIDENCE_TIERS)
  );
  const [pathways, setPathways] = useState<PathResult[]>([]);
  const [selectedIndex, setSelectedIndex] = useState(0);
  const [loading, setLoading] = useState(false);
  const [searched, setSearched] = useState(false);

  const graph = useMemo(() => buildGraph(reactions), []);

  const runSearch = useCallback(
    (src: string, tgt: string) => {
      if (!src || !tgt) return;
      setLoading(true);
      setSearched(true);
      // Use setTimeout to let the loading state render
      setTimeout(() => {
        const results = findKShortestPaths(graph, src, tgt, 10, {
          maxSteps,
          timeoutMs: 5000,
        });

        // Filter by evidence tiers
        const filtered = results.filter((p) => {
          return p.reactionIds.every((rId) => {
            const r = reactionMap.get(rId);
            return r ? enabledTiers.has(r.evidence_tier) : false;
          });
        });

        setPathways(filtered);
        setSelectedIndex(0);
        setLoading(false);
      }, 0);
    },
    [graph, maxSteps, enabledTiers]
  );

  // Run search when URL params change
  useEffect(() => {
    if (sourceId && targetId) {
      setSource(compoundMap.get(sourceId) ?? null);
      setTarget(compoundMap.get(targetId) ?? null);
      runSearch(sourceId, targetId);
    }
  }, [sourceId, targetId, runSearch]);

  function handleFind() {
    if (!source || !target) return;
    router.push(`/pathways?source=${source.id}&target=${target.id}`);
  }

  function toggleTier(tier: EvidenceTier) {
    setEnabledTiers((prev) => {
      const next = new Set(prev);
      if (next.has(tier)) {
        next.delete(tier);
      } else {
        next.add(tier);
      }
      return next;
    });
  }

  return (
    <div className="flex flex-1 flex-col px-4 py-6 sm:px-6">
      {/* Search bar */}
      <div className="mx-auto w-full max-w-4xl">
        <div className="flex items-end gap-3">
          <div className="flex-1">
            <label className="mb-1.5 block text-xs font-medium uppercase tracking-wider text-zinc-500">
              Source
            </label>
            <CompoundSearch
              placeholder="Starting compound..."
              value={source?.name ?? ""}
              onSelect={(c) => setSource(c)}
            />
          </div>
          <ArrowRight className="mb-2.5 h-5 w-5 shrink-0 text-zinc-600" />
          <div className="flex-1">
            <label className="mb-1.5 block text-xs font-medium uppercase tracking-wider text-zinc-500">
              Target
            </label>
            <CompoundSearch
              placeholder="Target compound..."
              value={target?.name ?? ""}
              onSelect={(c) => setTarget(c)}
            />
          </div>
          <button
            onClick={handleFind}
            disabled={!source || !target}
            className="mb-0.5 shrink-0 rounded-lg bg-zinc-100 px-4 py-2 text-sm font-medium text-zinc-900 transition-colors hover:bg-white disabled:cursor-not-allowed disabled:opacity-40"
          >
            Find
          </button>
        </div>

        {/* Filters */}
        <div className="mt-4 flex flex-wrap items-center gap-4">
          <div className="flex items-center gap-2">
            <label className="text-xs text-zinc-500">Max steps:</label>
            <input
              type="range"
              min={1}
              max={12}
              value={maxSteps}
              onChange={(e) => setMaxSteps(Number(e.target.value))}
              className="h-1.5 w-24 accent-zinc-400"
            />
            <span className="text-xs font-medium text-zinc-300">
              {maxSteps}
            </span>
          </div>

          <div className="flex items-center gap-2">
            <span className="text-xs text-zinc-500">Evidence:</span>
            {EVIDENCE_TIERS.map((tier) => (
              <label
                key={tier}
                className="flex items-center gap-1 text-xs text-zinc-400"
              >
                <input
                  type="checkbox"
                  checked={enabledTiers.has(tier)}
                  onChange={() => toggleTier(tier)}
                  className="h-3 w-3 accent-zinc-400"
                />
                <span className="capitalize">{tier}</span>
              </label>
            ))}
          </div>
        </div>
      </div>

      {/* Results */}
      <div className="mx-auto mt-6 w-full max-w-4xl flex-1">
        {loading ? (
          <div className="flex h-48 items-center justify-center">
            <Loader2 className="h-6 w-6 animate-spin text-zinc-500" />
            <span className="ml-2 text-sm text-zinc-500">
              Finding pathways...
            </span>
          </div>
        ) : !searched ? (
          <div className="flex h-48 items-center justify-center text-sm text-zinc-500">
            Enter source and target compounds to find pathways
          </div>
        ) : pathways.length === 0 ? (
          <div className="flex h-48 items-center justify-center text-sm text-zinc-500">
            No pathways found between these compounds
          </div>
        ) : (
          <div className="grid gap-6 lg:grid-cols-[45%_1fr]">
            <div>
              <PathwayList
                pathways={pathways}
                reactionMap={reactionMap}
                selectedIndex={selectedIndex}
                onSelect={setSelectedIndex}
              />
            </div>
            <div>
              {pathways[selectedIndex] && (
                <PathwayDetail
                  pathway={pathways[selectedIndex]}
                  reactionMap={reactionMap}
                />
              )}
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

export default function PathwayFinderPage() {
  return (
    <Suspense
      fallback={
        <div className="flex h-48 items-center justify-center">
          <Loader2 className="h-6 w-6 animate-spin text-zinc-500" />
        </div>
      }
    >
      <PathwayFinderContent />
    </Suspense>
  );
}
