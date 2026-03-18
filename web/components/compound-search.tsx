"use client";

import { useState, useRef, useEffect, useCallback } from "react";
import { searchCompounds } from "@/lib/search";
import { compounds } from "@/lib/data";
import { COMPOUND_TYPE_COLORS } from "@/lib/utils";
import type { Compound } from "@/lib/types";
import { Search } from "lucide-react";

interface CompoundSearchProps {
  onSelect: (compound: Compound) => void;
  placeholder?: string;
  value?: string;
}

export function CompoundSearch({
  onSelect,
  placeholder = "Search compounds...",
  value = "",
}: CompoundSearchProps) {
  const [query, setQuery] = useState(value);
  const [results, setResults] = useState<Compound[]>([]);
  const [open, setOpen] = useState(false);
  const [highlightIndex, setHighlightIndex] = useState(0);
  const containerRef = useRef<HTMLDivElement>(null);
  const inputRef = useRef<HTMLInputElement>(null);

  useEffect(() => {
    setQuery(value);
  }, [value]);

  const handleSearch = useCallback((q: string) => {
    setQuery(q);
    if (q.length < 1) {
      setResults([]);
      setOpen(false);
      return;
    }
    const hits = searchCompounds(q, compounds, 8);
    setResults(hits);
    setOpen(hits.length > 0);
    setHighlightIndex(0);
  }, []);

  useEffect(() => {
    function handleClickOutside(e: MouseEvent) {
      if (
        containerRef.current &&
        !containerRef.current.contains(e.target as Node)
      ) {
        setOpen(false);
      }
    }
    document.addEventListener("mousedown", handleClickOutside);
    return () => document.removeEventListener("mousedown", handleClickOutside);
  }, []);

  function select(compound: Compound) {
    setQuery(compound.name);
    setOpen(false);
    onSelect(compound);
  }

  function handleKeyDown(e: React.KeyboardEvent) {
    if (!open) return;
    if (e.key === "ArrowDown") {
      e.preventDefault();
      setHighlightIndex((i) => Math.min(i + 1, results.length - 1));
    } else if (e.key === "ArrowUp") {
      e.preventDefault();
      setHighlightIndex((i) => Math.max(i - 1, 0));
    } else if (e.key === "Enter") {
      e.preventDefault();
      if (results[highlightIndex]) {
        select(results[highlightIndex]);
      }
    } else if (e.key === "Escape") {
      setOpen(false);
    }
  }

  return (
    <div ref={containerRef} className="relative w-full">
      <div className="flex items-center gap-2 rounded-lg border border-zinc-800 bg-zinc-900/80 px-3 py-2 transition-colors focus-within:border-zinc-600">
        <Search className="h-4 w-4 shrink-0 text-zinc-500" />
        <input
          ref={inputRef}
          type="text"
          value={query}
          onChange={(e) => handleSearch(e.target.value)}
          onFocus={() => {
            if (query.length >= 1 && results.length > 0) setOpen(true);
          }}
          onKeyDown={handleKeyDown}
          placeholder={placeholder}
          className="w-full bg-transparent text-sm text-zinc-100 placeholder:text-zinc-500 outline-none"
        />
      </div>
      {open && results.length > 0 && (
        <div className="absolute z-50 mt-1 w-full rounded-lg border border-zinc-800 bg-zinc-900 shadow-xl">
          {results.map((compound, i) => (
            <button
              key={compound.id}
              type="button"
              className={`flex w-full items-center gap-3 px-3 py-2 text-left text-sm transition-colors ${
                i === highlightIndex
                  ? "bg-zinc-800 text-zinc-100"
                  : "text-zinc-300 hover:bg-zinc-800/50"
              } ${i === 0 ? "rounded-t-lg" : ""} ${
                i === results.length - 1 ? "rounded-b-lg" : ""
              }`}
              onMouseEnter={() => setHighlightIndex(i)}
              onClick={() => select(compound)}
            >
              <div className="flex-1">
                <span className="font-medium text-zinc-100">
                  {compound.name}
                </span>
                <span className="ml-2 text-xs text-zinc-500">
                  {compound.id}
                </span>
              </div>
              <span
                className={`text-xs ${COMPOUND_TYPE_COLORS[compound.type]}`}
              >
                {compound.type}
              </span>
              {compound.chirality !== "achiral" && (
                <span className="text-xs text-zinc-500">
                  {compound.chirality}
                </span>
              )}
            </button>
          ))}
        </div>
      )}
    </div>
  );
}
