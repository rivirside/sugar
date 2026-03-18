import Link from "next/link";
import {
  compounds,
  compoundMap,
  reactionMap,
  getReactionsForCompound,
} from "@/lib/data";
import { EvidenceBadge } from "@/components/evidence-badge";
import { COMPOUND_TYPE_COLORS, formatYield } from "@/lib/utils";
import {
  Table,
  TableHeader,
  TableBody,
  TableHead,
  TableRow,
  TableCell,
} from "@/components/ui/table";
import { ArrowRight } from "lucide-react";

export async function generateStaticParams() {
  return compounds.map((c) => ({ id: c.id }));
}

export default async function CompoundDetailPage({
  params,
}: {
  params: Promise<{ id: string }>;
}) {
  const { id } = await params;
  const compound = compoundMap.get(id);

  if (!compound) {
    return (
      <div className="flex h-64 items-center justify-center text-zinc-500">
        Compound not found: {id}
      </div>
    );
  }

  const relatedReactions = getReactionsForCompound(id);
  const incoming = relatedReactions.filter((r) => r.products.includes(id));
  const outgoing = relatedReactions.filter((r) => r.substrates.includes(id));

  return (
    <div className="flex flex-1 flex-col px-4 py-6 sm:px-6">
      <div className="mx-auto w-full max-w-4xl">
        {/* Breadcrumb */}
        <div className="flex items-center gap-2 text-xs text-zinc-500">
          <Link href="/compounds" className="hover:text-zinc-300">
            Compounds
          </Link>
          <span>/</span>
          <span className="text-zinc-300">{compound.name}</span>
        </div>

        {/* Header */}
        <div className="mt-4">
          <h1 className="text-3xl font-bold text-zinc-100">{compound.name}</h1>
          <div className="mt-2 flex flex-wrap items-center gap-3 text-sm">
            <span className="font-mono text-zinc-400">{compound.id}</span>
            <span
              className={`capitalize ${COMPOUND_TYPE_COLORS[compound.type]}`}
            >
              {compound.type.replace("_", " ")}
            </span>
            <span className="text-zinc-400">{compound.chirality}</span>
            <span className="font-mono text-zinc-400">{compound.formula}</span>
          </div>
        </div>

        {/* Details grid */}
        <div className="mt-6 grid gap-4 sm:grid-cols-2">
          <div className="rounded-lg border border-zinc-800 bg-zinc-900/60 p-4">
            <h3 className="text-xs font-medium uppercase tracking-wider text-zinc-500">
              Properties
            </h3>
            <dl className="mt-3 space-y-2 text-sm">
              <div className="flex justify-between">
                <dt className="text-zinc-500">Carbon count</dt>
                <dd className="text-zinc-200">{compound.carbons}</dd>
              </div>
              <div className="flex justify-between">
                <dt className="text-zinc-500">Chirality</dt>
                <dd className="text-zinc-200">{compound.chirality}</dd>
              </div>
              <div className="flex justify-between">
                <dt className="text-zinc-500">Commercial</dt>
                <dd className="text-zinc-200">
                  {compound.commercial ? "Yes" : "No"}
                </dd>
              </div>
              {compound.cost_usd_per_kg !== null && (
                <div className="flex justify-between">
                  <dt className="text-zinc-500">Cost (USD/kg)</dt>
                  <dd className="text-zinc-200">
                    ${compound.cost_usd_per_kg}
                  </dd>
                </div>
              )}
            </dl>
          </div>

          <div className="rounded-lg border border-zinc-800 bg-zinc-900/60 p-4">
            <h3 className="text-xs font-medium uppercase tracking-wider text-zinc-500">
              Stereocenters
            </h3>
            {compound.stereocenters.length > 0 ? (
              <div className="mt-3 flex flex-wrap gap-2">
                {compound.stereocenters.map((sc, i) => (
                  <span
                    key={i}
                    className="rounded-full bg-zinc-800 px-2.5 py-1 text-xs font-medium text-zinc-300"
                  >
                    C{i + 1}: {sc}
                  </span>
                ))}
              </div>
            ) : (
              <p className="mt-3 text-sm text-zinc-500">
                No stereocenters (achiral)
              </p>
            )}
          </div>
        </div>

        {/* Quick action */}
        <div className="mt-6">
          <Link
            href={`/pathways?source=${compound.id}`}
            className="inline-flex items-center gap-2 rounded-lg border border-zinc-800 bg-zinc-900/60 px-4 py-2 text-sm text-zinc-300 transition-colors hover:border-zinc-600 hover:text-zinc-100"
          >
            Find pathways from this compound
            <ArrowRight className="h-3.5 w-3.5" />
          </Link>
        </div>

        {/* Outgoing reactions */}
        {outgoing.length > 0 && (
          <div className="mt-8">
            <h2 className="text-lg font-semibold text-zinc-200">
              Outgoing Reactions ({outgoing.length})
            </h2>
            <div className="mt-3">
              <Table>
                <TableHeader>
                  <TableRow>
                    <TableHead>Product</TableHead>
                    <TableHead>Type</TableHead>
                    <TableHead>Evidence</TableHead>
                    <TableHead>Yield</TableHead>
                  </TableRow>
                </TableHeader>
                <TableBody>
                  {outgoing.map((r) => {
                    const product = compoundMap.get(r.products[0]);
                    return (
                      <TableRow key={r.id}>
                        <TableCell>
                          <Link
                            href={`/compound/${r.products[0]}`}
                            className="text-zinc-200 hover:text-white"
                          >
                            {product?.name ?? r.products[0]}
                          </Link>
                        </TableCell>
                        <TableCell className="text-xs capitalize text-zinc-400">
                          {r.reaction_type.replace("_", " ")}
                        </TableCell>
                        <TableCell>
                          <EvidenceBadge tier={r.evidence_tier} />
                        </TableCell>
                        <TableCell className="text-zinc-400">
                          {formatYield(r.yield)}
                        </TableCell>
                      </TableRow>
                    );
                  })}
                </TableBody>
              </Table>
            </div>
          </div>
        )}

        {/* Incoming reactions */}
        {incoming.length > 0 && (
          <div className="mt-8">
            <h2 className="text-lg font-semibold text-zinc-200">
              Incoming Reactions ({incoming.length})
            </h2>
            <div className="mt-3">
              <Table>
                <TableHeader>
                  <TableRow>
                    <TableHead>Substrate</TableHead>
                    <TableHead>Type</TableHead>
                    <TableHead>Evidence</TableHead>
                    <TableHead>Yield</TableHead>
                  </TableRow>
                </TableHeader>
                <TableBody>
                  {incoming.map((r) => {
                    const substrate = compoundMap.get(r.substrates[0]);
                    return (
                      <TableRow key={r.id}>
                        <TableCell>
                          <Link
                            href={`/compound/${r.substrates[0]}`}
                            className="text-zinc-200 hover:text-white"
                          >
                            {substrate?.name ?? r.substrates[0]}
                          </Link>
                        </TableCell>
                        <TableCell className="text-xs capitalize text-zinc-400">
                          {r.reaction_type.replace("_", " ")}
                        </TableCell>
                        <TableCell>
                          <EvidenceBadge tier={r.evidence_tier} />
                        </TableCell>
                        <TableCell className="text-zinc-400">
                          {formatYield(r.yield)}
                        </TableCell>
                      </TableRow>
                    );
                  })}
                </TableBody>
              </Table>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
