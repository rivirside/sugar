import { EVIDENCE_COLORS } from "@/lib/utils";
import type { EvidenceTier } from "@/lib/types";

interface EvidenceBadgeProps {
  tier: EvidenceTier;
  className?: string;
}

export function EvidenceBadge({ tier, className = "" }: EvidenceBadgeProps) {
  return (
    <span
      className={`inline-flex items-center rounded-full border px-2 py-0.5 text-xs font-medium capitalize ${EVIDENCE_COLORS[tier]} ${className}`}
    >
      {tier}
    </span>
  );
}
