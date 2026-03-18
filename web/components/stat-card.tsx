import Link from "next/link";

interface StatCardProps {
  label: string;
  value: string | number;
  href?: string;
  color?: string;
}

export function StatCard({
  label,
  value,
  href,
  color = "text-zinc-100",
}: StatCardProps) {
  const inner = (
    <div className="rounded-lg border border-zinc-800 bg-zinc-900/60 px-5 py-4 transition-colors hover:border-zinc-700">
      <p className="text-xs font-medium uppercase tracking-wider text-zinc-500">
        {label}
      </p>
      <p className={`mt-1 text-2xl font-semibold tabular-nums ${color}`}>
        {value}
      </p>
    </div>
  );

  if (href) {
    return (
      <Link href={href} className="block">
        {inner}
      </Link>
    );
  }

  return inner;
}
