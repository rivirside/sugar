"""ChEBI database importer.

Primary: bulk TSV download from ChEBI FTP.
Fallback: REST API per-compound.
Fallback on fallback: log warning and skip.
"""

import csv
import io
import logging
import os
import gzip

import requests

from pipeline.import_.cache import read_cache, write_cache, is_cache_fresh, write_raw_cache, read_raw_cache

logger = logging.getLogger(__name__)

CHEBI_FTP_COMPOUNDS = "https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/compounds.tsv.gz"
CHEBI_FTP_NAMES = "https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/names.tsv.gz"
CHEBI_REST_BASE = "https://www.ebi.ac.uk/webservices/chebi/2.0"


def fetch_chebi_bulk(cache_dir: str, refresh: bool = False) -> dict:
    """Download and parse ChEBI bulk TSV files. Returns a chebi_index dict."""
    index_cache = read_cache(cache_dir, "chebi", "index.json")
    if index_cache and not refresh and is_cache_fresh(cache_dir, "chebi", "index.json"):
        logger.info("Using cached ChEBI index (%d entries)", len(index_cache))
        return index_cache

    try:
        logger.info("Downloading ChEBI compounds TSV...")
        compounds_data = _download_tsv_gz(CHEBI_FTP_COMPOUNDS)
        logger.info("Downloading ChEBI names TSV...")
        names_data = _download_tsv_gz(CHEBI_FTP_NAMES)

        compounds = parse_chebi_compounds_tsv(compounds_data)
        names = parse_chebi_names_tsv(names_data)
        index = build_chebi_index(compounds, names)

        write_cache(cache_dir, "chebi", "index.json", index)
        logger.info("Built ChEBI index with %d entries", len(index))
        return index
    except Exception as e:
        logger.warning("ChEBI bulk download failed: %s. Falling back to cached data.", e)
        if index_cache:
            return index_cache
        return {}


def fetch_chebi_rest(compound_name: str) -> dict | None:
    """Fetch a single compound from ChEBI REST API. Returns entry dict or None."""
    try:
        search_url = f"{CHEBI_REST_BASE}/getLiteEntity?search={compound_name}&searchCategory=ALL&maximumResults=5"
        resp = requests.get(search_url, timeout=30)
        if resp.status_code == 200:
            return None  # XML parsing not implemented yet
        return None
    except Exception as e:
        logger.warning("ChEBI REST lookup failed for %s: %s", compound_name, e)
        return None


def parse_chebi_compounds_tsv(tsv_content: str) -> dict:
    """Parse ChEBI compounds.tsv. Returns {chebi_numeric_id: {"name": str, "chebi_id": str}}."""
    entries = {}
    reader = csv.DictReader(io.StringIO(tsv_content), delimiter="\t")
    for row in reader:
        if row.get("STATUS", "") != "C":
            continue
        chebi_num_id = row.get("ID", "").strip()
        name = row.get("NAME", "").strip()
        if chebi_num_id and name:
            entries[chebi_num_id] = {"name": name, "chebi_id": f"CHEBI:{chebi_num_id}"}
    return entries


def parse_chebi_names_tsv(tsv_content: str) -> dict:
    """Parse ChEBI names.tsv. Returns {chebi_numeric_id: [synonym1, synonym2, ...]}."""
    names: dict[str, list[str]] = {}
    reader = csv.DictReader(io.StringIO(tsv_content), delimiter="\t")
    for row in reader:
        compound_id = row.get("COMPOUND_ID", "").strip()
        name = row.get("NAME", "").strip()
        if compound_id and name:
            names.setdefault(compound_id, []).append(name)
    return names


def build_chebi_index(compounds: dict, names: dict) -> dict:
    """Build lookup index keyed by lowercase name/synonym -> ChEBI entry."""
    index = {}
    for chebi_num_id, compound_info in compounds.items():
        entry = {
            "chebi_id": compound_info["chebi_id"],
            "name": compound_info["name"],
            "synonyms": names.get(chebi_num_id, []),
            "formula": None,
            "inchi": None,
            "smiles": None,
            "kegg_id": None,
            "pubchem_id": None,
        }
        index[compound_info["name"].lower()] = entry
        for synonym in entry["synonyms"]:
            index[synonym.lower()] = entry
    return index


def _download_tsv_gz(url: str) -> str:
    """Download a gzipped TSV file and return decompressed content as string."""
    resp = requests.get(url, timeout=300, stream=True)
    resp.raise_for_status()
    content = gzip.decompress(resp.content)
    return content.decode("utf-8")
