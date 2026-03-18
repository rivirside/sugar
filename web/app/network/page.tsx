"use client";

import { useState, useCallback, useMemo, useRef } from "react";
import dynamic from "next/dynamic";
import { compounds, reactions, compoundMap, reactionMap } from "@/lib/data";
import { EvidenceBadge } from "@/components/evidence-badge";
import {
  COMPOUND_TYPE_COLORS,
  EVIDENCE_DOT_COLORS,
  formatYield,
} from "@/lib/utils";
import type { CompoundType, EvidenceTier } from "@/lib/types";
import type cytoscape from "cytoscape";
import Link from "next/link";
import { X } from "lucide-react";

const CytoscapeComponent = dynamic(() => import("react-cytoscapejs"), {
  ssr: false,
  loading: () => (
    <div className="flex h-full items-center justify-center text-sm text-zinc-500">
      Loading graph...
    </div>
  ),
});

const COMPOUND_TYPES: CompoundType[] = [
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

// Map compound types to hex colors for Cytoscape
const TYPE_HEX: Record<CompoundType, string> = {
  aldose: "#c084fc", // purple-400
  ketose: "#fb923c", // orange-400
  polyol: "#818cf8", // indigo-400
  phosphate: "#facc15", // yellow-400
  acid: "#f87171", // red-400
  lactone: "#f472b6", // pink-400
  amino_sugar: "#2dd4bf", // teal-400
  nucleotide_sugar: "#22d3ee", // cyan-400
  deoxy_sugar: "#a3e635", // lime-400
  disaccharide: "#fb7185", // rose-400
};

const EVIDENCE_HEX: Record<EvidenceTier, string> = {
  validated: "#4ade80",
  predicted: "#60a5fa",
  inferred: "#fbbf24",
  hypothetical: "#71717a",
};

interface InfoPanelData {
  type: "node" | "edge";
  id: string;
}

export default function NetworkPage() {
  const [enabledTypes, setEnabledTypes] = useState<Set<CompoundType>>(
    new Set(COMPOUND_TYPES)
  );
  const [carbonRange, setCarbonRange] = useState<[number, number]>([2, 12]);
  const [infoPanelData, setInfoPanelData] = useState<InfoPanelData | null>(
    null
  );
  const cyRef = useRef<cytoscape.Core | null>(null);

  function toggleType(type: CompoundType) {
    setEnabledTypes((prev) => {
      const next = new Set(prev);
      if (next.has(type)) next.delete(type);
      else next.add(type);
      return next;
    });
  }

  // Build Cytoscape elements from filtered data
  const elements = useMemo(() => {
    const filteredCompounds = compounds.filter(
      (c) =>
        enabledTypes.has(c.type) &&
        c.carbons >= carbonRange[0] &&
        c.carbons <= carbonRange[1]
    );
    const compoundIds = new Set(filteredCompounds.map((c) => c.id));

    const nodes: cytoscape.ElementDefinition[] = filteredCompounds.map((c) => ({
      data: {
        id: c.id,
        label: c.name,
        type: c.type,
        color: TYPE_HEX[c.type],
      },
    }));

    const edges: cytoscape.ElementDefinition[] = reactions
      .filter(
        (r) =>
          r.substrates.length === 1 &&
          r.products.length === 1 &&
          compoundIds.has(r.substrates[0]) &&
          compoundIds.has(r.products[0])
      )
      .map((r) => ({
        data: {
          id: r.id,
          source: r.substrates[0],
          target: r.products[0],
          color: EVIDENCE_HEX[r.evidence_tier],
          evidence: r.evidence_tier,
        },
      }));

    return [...nodes, ...edges];
  }, [enabledTypes, carbonRange]);

  const stylesheet: cytoscape.StylesheetStyle[] = useMemo(
    () => [
      {
        selector: "node",
        style: {
          label: "data(label)",
          "font-size": "8px",
          color: "#d4d4d8",
          "text-valign": "bottom",
          "text-margin-y": 4,
          "background-color": "data(color)",
          width: 16,
          height: 16,
        } as cytoscape.Css.Node,
      },
      {
        selector: "edge",
        style: {
          width: 1,
          "line-color": "data(color)",
          "target-arrow-color": "data(color)",
          "target-arrow-shape": "triangle",
          "arrow-scale": 0.6,
          "curve-style": "bezier",
          opacity: 0.6,
        } as cytoscape.Css.Edge,
      },
      {
        selector: "node:selected",
        style: {
          "border-width": 2,
          "border-color": "#ffffff",
        } as cytoscape.Css.Node,
      },
    ],
    []
  );

  const layout = useMemo(
    () => ({
      name: "cose" as const,
      animate: false,
      nodeDimensionsIncludeLabels: true,
      idealEdgeLength: 60,
      nodeRepulsion: 8000,
      gravity: 0.5,
    }),
    []
  );

  const handleCy = useCallback((cy: cytoscape.Core) => {
    cyRef.current = cy;
    cy.on("tap", "node", (evt) => {
      const nodeId = evt.target.id();
      setInfoPanelData({ type: "node", id: nodeId });
    });
    cy.on("tap", "edge", (evt) => {
      const edgeId = evt.target.id();
      setInfoPanelData({ type: "edge", id: edgeId });
    });
    cy.on("tap", (evt) => {
      if (evt.target === cy) {
        setInfoPanelData(null);
      }
    });
  }, []);

  // Info panel content
  const infoPanelContent = useMemo(() => {
    if (!infoPanelData) return null;

    if (infoPanelData.type === "node") {
      const compound = compoundMap.get(infoPanelData.id);
      if (!compound) return null;
      return (
        <div className="space-y-3">
          <div className="flex items-center justify-between">
            <h3 className="text-sm font-semibold text-zinc-100">
              {compound.name}
            </h3>
            <button
              onClick={() => setInfoPanelData(null)}
              className="text-zinc-500 hover:text-zinc-300"
            >
              <X className="h-4 w-4" />
            </button>
          </div>
          <div className="space-y-1 text-xs">
            <p className="text-zinc-500">
              ID: <span className="font-mono text-zinc-300">{compound.id}</span>
            </p>
            <p className="text-zinc-500">
              Type:{" "}
              <span
                className={`capitalize ${COMPOUND_TYPE_COLORS[compound.type]}`}
              >
                {compound.type.replace("_", " ")}
              </span>
            </p>
            <p className="text-zinc-500">
              Carbons:{" "}
              <span className="text-zinc-300">{compound.carbons}</span>
            </p>
            <p className="text-zinc-500">
              Chirality:{" "}
              <span className="text-zinc-300">{compound.chirality}</span>
            </p>
            <p className="text-zinc-500">
              Formula:{" "}
              <span className="font-mono text-zinc-300">
                {compound.formula}
              </span>
            </p>
          </div>
          <Link
            href={`/compound/${compound.id}`}
            className="inline-block text-xs text-blue-400 hover:text-blue-300"
          >
            View detail page
          </Link>
        </div>
      );
    } else {
      const reaction = reactionMap.get(infoPanelData.id);
      if (!reaction) return null;
      const substrate = compoundMap.get(reaction.substrates[0]);
      const product = compoundMap.get(reaction.products[0]);
      return (
        <div className="space-y-3">
          <div className="flex items-center justify-between">
            <h3 className="text-sm font-semibold text-zinc-100">Reaction</h3>
            <button
              onClick={() => setInfoPanelData(null)}
              className="text-zinc-500 hover:text-zinc-300"
            >
              <X className="h-4 w-4" />
            </button>
          </div>
          <div className="space-y-1 text-xs">
            <p className="text-zinc-300">
              {substrate?.name ?? reaction.substrates[0]} {"\u2192"}{" "}
              {product?.name ?? reaction.products[0]}
            </p>
            <p className="text-zinc-500">
              Type:{" "}
              <span className="capitalize text-zinc-300">
                {reaction.reaction_type.replace("_", " ")}
              </span>
            </p>
            <div className="flex items-center gap-1 text-zinc-500">
              Evidence: <EvidenceBadge tier={reaction.evidence_tier} />
            </div>
            <p className="text-zinc-500">
              Yield:{" "}
              <span className="text-zinc-300">
                {formatYield(reaction.yield)}
              </span>
            </p>
            <p className="text-zinc-500">
              Cost:{" "}
              <span className="text-zinc-300">
                {reaction.cost_score.toFixed(2)}
              </span>
            </p>
          </div>
          <Link
            href={`/reaction/${reaction.id}`}
            className="inline-block text-xs text-blue-400 hover:text-blue-300"
          >
            View detail page
          </Link>
        </div>
      );
    }
  }, [infoPanelData]);

  return (
    <div className="flex flex-1 flex-col">
      {/* Mobile message */}
      <div className="block p-6 text-center text-sm text-zinc-500 sm:hidden">
        This network visualization is best viewed on desktop.
      </div>

      <div className="hidden flex-1 sm:flex">
        {/* Filter panel */}
        <div className="w-56 shrink-0 space-y-4 border-r border-zinc-800 p-4">
          <h2 className="text-sm font-semibold text-zinc-200">Filters</h2>

          {/* Compound type checkboxes */}
          <div className="space-y-2">
            <h3 className="text-xs font-medium uppercase tracking-wider text-zinc-500">
              Compound Type
            </h3>
            {COMPOUND_TYPES.map((type) => (
              <label
                key={type}
                className="flex items-center gap-2 text-xs text-zinc-400"
              >
                <input
                  type="checkbox"
                  checked={enabledTypes.has(type)}
                  onChange={() => toggleType(type)}
                  className="h-3 w-3 accent-zinc-400"
                />
                <span
                  className="h-2 w-2 rounded-full"
                  style={{ backgroundColor: TYPE_HEX[type] }}
                />
                <span className="capitalize">{type.replace("_", " ")}</span>
              </label>
            ))}
          </div>

          {/* Carbon range slider */}
          <div className="space-y-2">
            <h3 className="text-xs font-medium uppercase tracking-wider text-zinc-500">
              Carbon Range
            </h3>
            <div className="flex items-center gap-2">
              <input
                type="range"
                min={2}
                max={12}
                value={carbonRange[0]}
                onChange={(e) =>
                  setCarbonRange([Number(e.target.value), carbonRange[1]])
                }
                className="h-1 w-full accent-zinc-400"
              />
            </div>
            <div className="flex items-center gap-2">
              <input
                type="range"
                min={2}
                max={12}
                value={carbonRange[1]}
                onChange={(e) =>
                  setCarbonRange([carbonRange[0], Number(e.target.value)])
                }
                className="h-1 w-full accent-zinc-400"
              />
            </div>
            <p className="text-xs text-zinc-400">
              {carbonRange[0]} &ndash; {carbonRange[1]} carbons
            </p>
          </div>

          {/* Legend */}
          <div className="space-y-2">
            <h3 className="text-xs font-medium uppercase tracking-wider text-zinc-500">
              Edge Legend
            </h3>
            {(
              ["validated", "predicted", "inferred", "hypothetical"] as const
            ).map((tier) => (
              <div key={tier} className="flex items-center gap-2 text-xs">
                <span
                  className={`h-2 w-2 rounded-full ${EVIDENCE_DOT_COLORS[tier]}`}
                />
                <span className="capitalize text-zinc-400">{tier}</span>
              </div>
            ))}
          </div>

          <div className="pt-2 text-xs text-zinc-600">
            {elements.filter((e) => !e.data.source).length} nodes,{" "}
            {elements.filter((e) => e.data.source).length} edges
          </div>
        </div>

        {/* Graph area */}
        <div className="relative flex-1">
          <CytoscapeComponent
            elements={elements}
            stylesheet={stylesheet}
            layout={layout}
            style={{ width: "100%", height: "100%" }}
            cy={handleCy}
            minZoom={0.2}
            maxZoom={3}
            userZoomingEnabled
            userPanningEnabled
            boxSelectionEnabled={false}
          />

          {/* Info panel */}
          {infoPanelContent && (
            <div className="absolute right-4 top-4 w-64 rounded-lg border border-zinc-800 bg-zinc-900/95 p-4 shadow-xl backdrop-blur">
              {infoPanelContent}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
