"""Tests for pipeline.import_.cache module."""

import json
import os
import tempfile
import time

import pytest

from pipeline.import_.cache import read_cache, write_cache, is_cache_fresh, clear_cache


@pytest.fixture
def cache_dir(tmp_path):
    return str(tmp_path / "cache")


def test_write_and_read_cache(cache_dir):
    data = {"chebi_id": "CHEBI:17634", "name": "D-Glucose"}
    write_cache(cache_dir, "chebi", "glucose.json", data)
    result = read_cache(cache_dir, "chebi", "glucose.json")
    assert result == data


def test_read_missing_cache(cache_dir):
    result = read_cache(cache_dir, "chebi", "nonexistent.json")
    assert result is None


def test_is_cache_fresh_when_fresh(cache_dir):
    write_cache(cache_dir, "kegg", "C00031.json", {"id": "C00031"})
    assert is_cache_fresh(cache_dir, "kegg", "C00031.json", max_age_days=7)


def test_is_cache_fresh_when_missing(cache_dir):
    assert not is_cache_fresh(cache_dir, "kegg", "missing.json", max_age_days=7)


def test_clear_cache(cache_dir):
    write_cache(cache_dir, "chebi", "test.json", {"a": 1})
    clear_cache(cache_dir, "chebi")
    assert read_cache(cache_dir, "chebi", "test.json") is None


def test_write_cache_creates_subdirs(cache_dir):
    write_cache(cache_dir, "brenda", "1.2.3.4.json", {"ec": "1.2.3.4"})
    path = os.path.join(cache_dir, "brenda", "1.2.3.4.json")
    assert os.path.exists(path)
