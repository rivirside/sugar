"use client";

import { useState, useMemo } from "react";
import { useRouter } from "next/navigation";
import { compounds } from "@/lib/data";
import { COMPOUND_TYPE_COLORS } from "@/lib/utils";
import type { CompoundType, Chirality } from "@/lib/types";
import {
  Table,
  TableHeader,
  TableBody,
  TableHead,
  TableRow,
  TableCell,
} from "@/components/ui/table";
import { Search } from "lucide-react";

const ALL_TYPES: CompoundType[] = [
  "aldose",
  "ketose",
  "polyol",
  "phosphate",
  "acid",
  "lactone",
  "amino_sugar",
  "nucleotide_sugar",
  "deoxy_sugar",
  "disaccharide",
];

const ALL_CHIRALITIES: Chirality[] = ["D", "L", "achiral"];

export default function CompoundBrowserPage() {
  const router = useRouter();
  const [query, setQuery] = useState("");
  const [typeFilter, setTypeFilter] = useState<Set<CompoundType>>(new Set());
  const [chiralityFilter, setChiralityFilter] = useState<Set<Chirality>>(
    new Set()
  );
  const [carbonRange, setCarbonRange] = useState<[number, number]>([2, 12]);

  const filtered = useMemo(() => {
    return compounds.filter((c) => {
      // Text search
      if (query) {
        const q = query.toLowerCase();
        const matchesName = c.name.toLowerCase().includes(q);
        const matchesId = c.id.toLowerCase().includes(q);
        const matchesAlias = c.aliases.some((a) =>
          a.toLowerCase().includes(q)
        );
        if (!matchesName && !matchesId && !matchesAlias) return false;
      }

      // Type filter
      if (typeFilter.size > 0 && !typeFilter.has(c.type)) return false;

      // Chirality filter
      if (chiralityFilter.size > 0 && !chiralityFilter.has(c.chirality))
        return false;

      // Carbon range
      if (c.carbons < carbonRange[0] || c.carbons > carbonRange[1])
        return false;

      return true;
    });
  }, [query, typeFilter, chiralityFilter, carbonRange]);

  function toggleType(type: CompoundType) {
    setTypeFilter((prev) => {
      const next = new Set(prev);
      if (next.has(type)) {
        next.delete(type);
      } else {
        next.add(type);
      }
      return next;
    });
  }

  function toggleChirality(chirality: Chirality) {
    setChiralityFilter((prev) => {
      const next = new Set(prev);
      if (next.has(chirality)) {
        next.delete(chirality);
      } else {
        next.add(chirality);
      }
      return next;
    });
  }

  return (
    <div className="flex flex-1 flex-col px-4 py-6 sm:px-6">
      <div className="mx-auto w-full max-w-5xl">
        <h1 className="text-2xl font-bold text-zinc-100">Compounds</h1>
        <p className="mt-1 text-sm text-zinc-400">
          Browse {compounds.length} sugar metabolites
        </p>

        {/* Search + Filters */}
        <div className="mt-6 space-y-4">
          <div className="flex items-center gap-2 rounded-lg border border-zinc-800 bg-zinc-900/80 px-3 py-2">
            <Search className="h-4 w-4 text-zinc-500" />
            <input
              type="text"
              value={query}
              onChange={(e) => setQuery(e.target.value)}
              placeholder="Search by name, ID, or alias..."
              className="w-full bg-transparent text-sm text-zinc-100 placeholder:text-zinc-500 outline-none"
            />
          </div>

          <div className="flex flex-wrap gap-4">
            {/* Type filter */}
            <div className="flex flex-wrap items-center gap-2">
              <span className="text-xs text-zinc-500">Type:</span>
              {ALL_TYPES.map((type) => (
                <button
                  key={type}
                  onClick={() => toggleType(type)}
                  className={`rounded-full border px-2 py-0.5 text-xs capitalize transition-colors ${
                    typeFilter.has(type)
                      ? "border-zinc-600 bg-zinc-800 text-zinc-200"
                      : "border-zinc-800 text-zinc-500 hover:text-zinc-300"
                  }`}
                >
                  {type.replace("_", " ")}
                </button>
              ))}
            </div>

            {/* Chirality filter */}
            <div className="flex items-center gap-2">
              <span className="text-xs text-zinc-500">Chirality:</span>
              {ALL_CHIRALITIES.map((c) => (
                <button
                  key={c}
                  onClick={() => toggleChirality(c)}
                  className={`rounded-full border px-2 py-0.5 text-xs transition-colors ${
                    chiralityFilter.has(c)
                      ? "border-zinc-600 bg-zinc-800 text-zinc-200"
                      : "border-zinc-800 text-zinc-500 hover:text-zinc-300"
                  }`}
                >
                  {c}
                </button>
              ))}
            </div>

            {/* Carbon range */}
            <div className="flex items-center gap-2">
              <span className="text-xs text-zinc-500">Carbons:</span>
              <input
                type="range"
                min={2}
                max={12}
                value={carbonRange[0]}
                onChange={(e) =>
                  setCarbonRange([Number(e.target.value), carbonRange[1]])
                }
                className="h-1.5 w-16 accent-zinc-400"
              />
              <span className="text-xs text-zinc-300">
                {carbonRange[0]}&ndash;{carbonRange[1]}
              </span>
              <input
                type="range"
                min={2}
                max={12}
                value={carbonRange[1]}
                onChange={(e) =>
                  setCarbonRange([carbonRange[0], Number(e.target.value)])
                }
                className="h-1.5 w-16 accent-zinc-400"
              />
            </div>
          </div>
        </div>

        <p className="mt-4 text-xs text-zinc-500">{filtered.length} results</p>

        {/* Table */}
        <div className="mt-2">
          <Table>
            <TableHeader>
              <TableRow>
                <TableHead>Name</TableHead>
                <TableHead>Type</TableHead>
                <TableHead>Carbons</TableHead>
                <TableHead>Chirality</TableHead>
                <TableHead>Formula</TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {filtered.map((c) => (
                <TableRow
                  key={c.id}
                  className="cursor-pointer"
                  onClick={() => router.push(`/compound/${c.id}`)}
                >
                  <TableCell>
                    <span className="font-medium text-zinc-200">{c.name}</span>
                    <span className="ml-2 text-xs text-zinc-500">{c.id}</span>
                  </TableCell>
                  <TableCell>
                    <span
                      className={`text-xs capitalize ${COMPOUND_TYPE_COLORS[c.type]}`}
                    >
                      {c.type.replace("_", " ")}
                    </span>
                  </TableCell>
                  <TableCell className="text-zinc-300">{c.carbons}</TableCell>
                  <TableCell className="text-zinc-300">{c.chirality}</TableCell>
                  <TableCell className="font-mono text-xs text-zinc-400">
                    {c.formula}
                  </TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        </div>
      </div>
    </div>
  );
}
