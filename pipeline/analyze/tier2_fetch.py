"""Tier 2 enzyme index enrichment via UniProt and PDB APIs.

Fetches protein family size, UniProt accessions, and PDB structure counts
for each EC number in the enzyme index. Results are cached to avoid
redundant API calls on subsequent pipeline runs.
"""

import json
import time
from pathlib import Path

import requests

# Cache directory for Tier 2 data
CACHE_DIR = Path(__file__).resolve().parent.parent / "cache" / "tier2"

# Rate limiting: seconds between API calls
_RATE_LIMIT = 0.2


def _cache_path(ec_number: str) -> Path:
    """Return the cache file path for an EC number."""
    safe_name = ec_number.replace(".", "_")
    return CACHE_DIR / f"{ec_number}.json"


def _read_cache(ec_number: str) -> dict | None:
    """Read cached Tier 2 data for an EC number, if it exists."""
    path = _cache_path(ec_number)
    if path.exists():
        try:
            return json.loads(path.read_text())
        except (json.JSONDecodeError, OSError):
            return None
    return None


def _write_cache(ec_number: str, data: dict) -> None:
    """Write Tier 2 data to cache."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path = _cache_path(ec_number)
    path.write_text(json.dumps(data, indent=2))


def fetch_uniprot_ec(ec_number: str) -> dict:
    """Fetch reviewed UniProt entries for an EC number.

    Queries the UniProt REST API for Swiss-Prot (reviewed) entries
    annotated with this EC number.

    Args:
        ec_number: EC number string (e.g., "5.1.3.2").

    Returns:
        Dict with keys: family_size (int), uniprot_ids (list[str]),
        organisms (list[str]).
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"ec:{ec_number} AND reviewed:true",
        "format": "json",
        "fields": "accession,organism_name",
        "size": 500,
    }

    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        entries = data.get("results", [])

        uniprot_ids = []
        organisms = []
        for entry in entries:
            acc = entry.get("primaryAccession")
            if acc:
                uniprot_ids.append(acc)
            org = entry.get("organism", {}).get("scientificName")
            if org and org not in organisms:
                organisms.append(org)

        return {
            "family_size": len(uniprot_ids),
            "uniprot_ids": uniprot_ids,
            "organisms": organisms,
        }
    except Exception:
        return {
            "family_size": 0,
            "uniprot_ids": [],
            "organisms": [],
        }


def fetch_pdb_ec(ec_number: str) -> dict:
    """Fetch PDB structures for an EC number.

    Queries the RCSB PDB search API for structures annotated
    with this EC number.

    Args:
        ec_number: EC number string (e.g., "5.1.3.2").

    Returns:
        Dict with keys: pdb_count (int), pdb_ids (list[str]).
    """
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity.rcsb_ec_lineage.id",
                "operator": "exact_match",
                "value": ec_number,
            },
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "return_all_hits": True,
        },
    }

    try:
        resp = requests.post(url, json=query, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        result_set = data.get("result_set", [])
        pdb_ids = [r["identifier"] for r in result_set if "identifier" in r]

        return {
            "pdb_count": data.get("total_count", len(pdb_ids)),
            "pdb_ids": pdb_ids,
        }
    except Exception:
        return {
            "pdb_count": 0,
            "pdb_ids": [],
        }


def enrich_enzyme_index(
    index: dict,
    use_cache: bool = True,
) -> dict:
    """Enrich a Tier 1 enzyme index with UniProt and PDB data.

    For each EC number in the index, fetches family size, UniProt IDs,
    and PDB structure counts. Uses cached data when available.

    Args:
        index: Tier 1 enzyme index (keyed by EC number).
        use_cache: Whether to use/write the disk cache.

    Returns:
        Enriched copy of the index with Tier 2 fields populated.
    """
    enriched = {}

    for ec, entry in index.items():
        enriched_entry = dict(entry)

        # Check cache first
        cached = _read_cache(ec) if use_cache else None

        if cached:
            tier2 = cached
        else:
            # Fetch from APIs
            uniprot = fetch_uniprot_ec(ec)
            time.sleep(_RATE_LIMIT)
            pdb = fetch_pdb_ec(ec)
            time.sleep(_RATE_LIMIT)

            tier2 = {
                "family_size": uniprot["family_size"],
                "uniprot_ids": uniprot["uniprot_ids"],
                "organisms": uniprot["organisms"],
                "pdb_count": pdb["pdb_count"],
                "pdb_ids": pdb["pdb_ids"],
            }

            # Write to cache
            if use_cache and (tier2["family_size"] > 0 or tier2["pdb_count"] > 0):
                _write_cache(ec, tier2)

        # Merge Tier 2 data into entry
        if tier2["family_size"] > 0:
            enriched_entry["family_size"] = tier2["family_size"]
            enriched_entry["uniprot_ids"] = tier2["uniprot_ids"]
            # Merge organisms (preserve originals, add new)
            existing_orgs = set(enriched_entry.get("organisms", []))
            for org in tier2.get("organisms", []):
                existing_orgs.add(org)
            enriched_entry["organisms"] = sorted(existing_orgs)
        # If family_size is 0, keep original None (no data, not confirmed zero)

        if tier2["pdb_count"] > 0:
            enriched_entry["pdb_count"] = tier2["pdb_count"]
        # Keep None if no PDB data found

        enriched[ec] = enriched_entry

    return enriched
