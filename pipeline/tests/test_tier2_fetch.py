"""Tests for Tier 2 enzyme index enrichment (UniProt/PDB)."""

import json
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

from pipeline.analyze.tier2_fetch import (
    fetch_uniprot_ec,
    fetch_pdb_ec,
    enrich_enzyme_index,
    _cache_path,
)


# --- UniProt fetch tests ---

def _mock_uniprot_response(entries):
    """Build a mock requests.Response for UniProt REST API."""
    resp = MagicMock()
    resp.status_code = 200
    resp.json.return_value = {"results": entries}
    return resp


def test_fetch_uniprot_ec_parses_results():
    entries = [
        {
            "primaryAccession": "P12345",
            "organism": {"scientificName": "Escherichia coli"},
        },
        {
            "primaryAccession": "Q67890",
            "organism": {"scientificName": "Homo sapiens"},
        },
    ]
    with patch("pipeline.analyze.tier2_fetch.requests.get",
               return_value=_mock_uniprot_response(entries)):
        result = fetch_uniprot_ec("5.1.3.2")

    assert result["family_size"] == 2
    assert set(result["uniprot_ids"]) == {"P12345", "Q67890"}
    assert "Escherichia coli" in result["organisms"]
    assert "Homo sapiens" in result["organisms"]


def test_fetch_uniprot_ec_empty_results():
    with patch("pipeline.analyze.tier2_fetch.requests.get",
               return_value=_mock_uniprot_response([])):
        result = fetch_uniprot_ec("9.9.9.9")

    assert result["family_size"] == 0
    assert result["uniprot_ids"] == []
    assert result["organisms"] == []


def test_fetch_uniprot_ec_handles_network_error():
    with patch("pipeline.analyze.tier2_fetch.requests.get",
               side_effect=Exception("Connection refused")):
        result = fetch_uniprot_ec("5.1.3.2")

    assert result["family_size"] == 0
    assert result["uniprot_ids"] == []


# --- PDB fetch tests ---

def _mock_pdb_response(pdb_ids):
    """Build a mock requests.Response for RCSB PDB search API."""
    resp = MagicMock()
    resp.status_code = 200
    resp.json.return_value = {
        "total_count": len(pdb_ids),
        "result_set": [{"identifier": pid} for pid in pdb_ids],
    }
    return resp


def test_fetch_pdb_ec_parses_results():
    with patch("pipeline.analyze.tier2_fetch.requests.post",
               return_value=_mock_pdb_response(["1ABC", "2DEF", "3GHI"])):
        result = fetch_pdb_ec("5.1.3.2")

    assert result["pdb_count"] == 3
    assert set(result["pdb_ids"]) == {"1ABC", "2DEF", "3GHI"}


def test_fetch_pdb_ec_empty_results():
    resp = MagicMock()
    resp.status_code = 200
    resp.json.return_value = {"total_count": 0, "result_set": []}
    with patch("pipeline.analyze.tier2_fetch.requests.post", return_value=resp):
        result = fetch_pdb_ec("9.9.9.9")

    assert result["pdb_count"] == 0
    assert result["pdb_ids"] == []


def test_fetch_pdb_ec_handles_network_error():
    with patch("pipeline.analyze.tier2_fetch.requests.post",
               side_effect=Exception("Timeout")):
        result = fetch_pdb_ec("5.1.3.2")

    assert result["pdb_count"] == 0
    assert result["pdb_ids"] == []


# --- Enrichment integration tests ---

def test_enrich_enzyme_index_merges_tier2_data():
    index = {
        "5.1.3.2": {
            "name": "UDP-glucose 4-epimerase",
            "organisms": ["E. coli"],
            "known_substrates": ["D-GLC"],
            "reaction_count": 2,
            "family_size": None,
            "pdb_count": None,
            "uniprot_ids": None,
        }
    }

    uniprot_data = {
        "family_size": 47,
        "uniprot_ids": ["P12345", "Q67890"],
        "organisms": ["E. coli", "H. sapiens"],
    }
    pdb_data = {
        "pdb_count": 12,
        "pdb_ids": ["1ABC", "2DEF"],
    }

    with patch("pipeline.analyze.tier2_fetch.fetch_uniprot_ec",
               return_value=uniprot_data), \
         patch("pipeline.analyze.tier2_fetch.fetch_pdb_ec",
               return_value=pdb_data):
        enriched = enrich_enzyme_index(index, use_cache=False)

    entry = enriched["5.1.3.2"]
    assert entry["family_size"] == 47
    assert entry["pdb_count"] == 12
    assert "P12345" in entry["uniprot_ids"]
    # Original organisms preserved alongside new ones
    assert "E. coli" in entry["organisms"]
    assert "H. sapiens" in entry["organisms"]


def test_enrich_enzyme_index_uses_cache(tmp_path):
    index = {
        "5.1.3.2": {
            "name": "test",
            "organisms": [],
            "known_substrates": [],
            "reaction_count": 1,
            "family_size": None,
            "pdb_count": None,
            "uniprot_ids": None,
        }
    }

    # Write cache file
    cache_dir = tmp_path / "tier2"
    cache_dir.mkdir()
    cached = {
        "family_size": 10,
        "uniprot_ids": ["CACHED1"],
        "organisms": ["Cached org"],
        "pdb_count": 3,
        "pdb_ids": ["XABC"],
    }
    (cache_dir / "5.1.3.2.json").write_text(json.dumps(cached))

    with patch("pipeline.analyze.tier2_fetch.CACHE_DIR", cache_dir):
        enriched = enrich_enzyme_index(index, use_cache=True)

    entry = enriched["5.1.3.2"]
    assert entry["family_size"] == 10
    assert entry["pdb_count"] == 3
    assert "CACHED1" in entry["uniprot_ids"]


def test_enrich_enzyme_index_skips_failed_fetches():
    """If both API calls fail, Tier 1 data is preserved unchanged."""
    index = {
        "5.1.3.2": {
            "name": "test",
            "organisms": ["E. coli"],
            "known_substrates": ["D-GLC"],
            "reaction_count": 1,
            "family_size": None,
            "pdb_count": None,
            "uniprot_ids": None,
        }
    }

    fail_data = {"family_size": 0, "uniprot_ids": [], "organisms": []}
    pdb_fail = {"pdb_count": 0, "pdb_ids": []}

    with patch("pipeline.analyze.tier2_fetch.fetch_uniprot_ec",
               return_value=fail_data), \
         patch("pipeline.analyze.tier2_fetch.fetch_pdb_ec",
               return_value=pdb_fail):
        enriched = enrich_enzyme_index(index, use_cache=False)

    entry = enriched["5.1.3.2"]
    # family_size stays None when API returned 0 (no data, not "zero enzymes")
    assert entry["family_size"] is None
    assert entry["organisms"] == ["E. coli"]  # original preserved
