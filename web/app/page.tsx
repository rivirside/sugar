"use client";

import { useState } from "react";
import { useRouter } from "next/navigation";
import { CompoundSearch } from "@/components/compound-search";
import { StatCard } from "@/components/stat-card";
import { metadata as pipelineMetadata } from "@/lib/data";
import type { Compound } from "@/lib/types";
import { ArrowRight } from "lucide-react";

const POPULAR_SEARCHES = [
  { source: "D-GLC", target: "D-MAN", label: "Glucose \u2192 Mannose" },
  { source: "D-GLC", target: "D-GAL", label: "Glucose \u2192 Galactose" },
  { source: "D-FRU", target: "D-SOR", label: "Fructose \u2192 Sorbitol" },
  { source: "D-MAN", target: "L-GUL", label: "Mannose \u2192 L-Gulose" },
];

export default function DashboardPage() {
  const router = useRouter();
  const [source, setSource] = useState<Compound | null>(null);
  const [target, setTarget] = useState<Compound | null>(null);

  function handleFindPathways() {
    if (!source || !target) return;
    router.push(`/pathways?source=${source.id}&target=${target.id}`);
  }

  function handlePopular(src: string, tgt: string) {
    router.push(`/pathways?source=${src}&target=${tgt}`);
  }

  const counts = pipelineMetadata.counts;

  return (
    <div className="flex flex-1 flex-col items-center px-4 pt-16 pb-12">
      {/* Hero */}
      <div className="w-full max-w-2xl text-center">
        <h1 className="text-4xl font-bold tracking-tight text-zinc-100">
          Find a synthesis pathway
        </h1>
        <p className="mt-3 text-zinc-400">
          Explore enzymatic routes between sugar metabolites
        </p>
      </div>

      {/* Search inputs */}
      <div className="mt-10 flex w-full max-w-2xl items-center gap-3">
        <div className="flex-1">
          <label className="mb-1.5 block text-xs font-medium uppercase tracking-wider text-zinc-500">
            Source
          </label>
          <CompoundSearch
            placeholder="Starting compound..."
            onSelect={setSource}
          />
        </div>

        <ArrowRight className="mt-6 h-5 w-5 shrink-0 text-zinc-600" />

        <div className="flex-1">
          <label className="mb-1.5 block text-xs font-medium uppercase tracking-wider text-zinc-500">
            Target
          </label>
          <CompoundSearch
            placeholder="Target compound..."
            onSelect={setTarget}
          />
        </div>
      </div>

      {/* Find button */}
      <button
        onClick={handleFindPathways}
        disabled={!source || !target}
        className="mt-6 rounded-lg bg-zinc-100 px-6 py-2.5 text-sm font-medium text-zinc-900 transition-colors hover:bg-white disabled:cursor-not-allowed disabled:opacity-40"
      >
        Find Pathways
      </button>

      {/* Popular searches */}
      <div className="mt-8 flex flex-wrap justify-center gap-2">
        {POPULAR_SEARCHES.map((p) => (
          <button
            key={`${p.source}-${p.target}`}
            onClick={() => handlePopular(p.source, p.target)}
            className="rounded-full border border-zinc-800 bg-zinc-900/60 px-3 py-1.5 text-xs text-zinc-400 transition-colors hover:border-zinc-600 hover:text-zinc-200"
          >
            {p.label}
          </button>
        ))}
      </div>

      {/* Stats row */}
      <div className="mt-16 grid w-full max-w-3xl grid-cols-2 gap-4 sm:grid-cols-4">
        <StatCard
          label="Compounds"
          value={counts.total_compounds}
          href="/compounds"
          color="text-purple-400"
        />
        <StatCard
          label="Reactions"
          value={counts.total_reactions}
          href="/reactions"
          color="text-blue-400"
        />
        <StatCard
          label="Monosaccharides"
          value={counts.monosaccharides ?? 0}
          href="/compounds"
          color="text-orange-400"
        />
        <StatCard
          label="Polyols"
          value={counts.polyols ?? 0}
          href="/compounds"
          color="text-indigo-400"
        />
      </div>
    </div>
  );
}
