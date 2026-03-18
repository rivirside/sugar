"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";
import { Search } from "lucide-react";

const NAV_LINKS = [
  { href: "/", label: "Dashboard" },
  { href: "/pathways", label: "Pathways" },
  { href: "/compounds", label: "Compounds" },
  { href: "/reactions", label: "Reactions" },
  { href: "/network", label: "Network" },
  { href: "/about", label: "About" },
];

export function NavBar() {
  const pathname = usePathname();

  return (
    <nav className="sticky top-0 z-50 flex h-12 items-center gap-6 border-b border-[#1e1e1e] bg-[#0a0a0a]/90 px-6 backdrop-blur">
      {/* Wordmark */}
      <Link href="/" className="mr-2 text-sm font-semibold tracking-widest text-white">
        SUGAR
      </Link>

      {/* Nav links */}
      <div className="flex items-center gap-1">
        {NAV_LINKS.map(({ href, label }) => {
          const isActive = href === "/" ? pathname === "/" : pathname.startsWith(href);
          return (
            <Link
              key={href}
              href={href}
              className={[
                "rounded px-3 py-1.5 text-sm transition-colors",
                isActive
                  ? "bg-white/10 text-white"
                  : "text-zinc-400 hover:text-zinc-200",
              ].join(" ")}
            >
              {label}
            </Link>
          );
        })}
      </div>

      {/* Spacer */}
      <div className="flex-1" />

      {/* Search button */}
      <button
        className="flex items-center gap-2 rounded border border-[#1e1e1e] bg-[#141414] px-3 py-1.5 text-sm text-zinc-500 transition-colors hover:border-zinc-600 hover:text-zinc-300"
        aria-label="Open search"
      >
        <Search className="h-3.5 w-3.5" />
        <span>Search...</span>
        <kbd className="ml-1 rounded bg-[#1e1e1e] px-1.5 py-0.5 font-mono text-xs text-zinc-500">
          ⌘K
        </kbd>
      </button>
    </nav>
  );
}
