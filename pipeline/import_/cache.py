"""Shared caching utilities for import pipeline.

Provides read/write/freshness-check for cached API responses.
All cached data lives in pipeline/cache/ (gitignored).
"""

import json
import os
import shutil
import time


def _cache_path(cache_dir: str, source: str, filename: str) -> str:
    return os.path.join(cache_dir, source, filename)


def write_cache(cache_dir: str, source: str, filename: str, data: dict | list) -> str:
    """Write data to cache as JSON. Returns the file path."""
    path = _cache_path(cache_dir, source, filename)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    return path


def read_cache(cache_dir: str, source: str, filename: str) -> dict | list | None:
    """Read cached data. Returns None if file doesn't exist."""
    path = _cache_path(cache_dir, source, filename)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def is_cache_fresh(cache_dir: str, source: str, filename: str, max_age_days: int = 30) -> bool:
    """Check if a cached file exists and is younger than max_age_days."""
    path = _cache_path(cache_dir, source, filename)
    if not os.path.exists(path):
        return False
    mtime = os.path.getmtime(path)
    age_days = (time.time() - mtime) / 86400
    return age_days < max_age_days


def clear_cache(cache_dir: str, source: str) -> None:
    """Remove all cached files for a specific source."""
    source_dir = os.path.join(cache_dir, source)
    if os.path.exists(source_dir):
        shutil.rmtree(source_dir)


def write_raw_cache(cache_dir: str, source: str, filename: str, content: bytes) -> str:
    """Write raw bytes to cache (for TSV/binary downloads). Returns file path."""
    path = _cache_path(cache_dir, source, filename)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        f.write(content)
    return path


def read_raw_cache(cache_dir: str, source: str, filename: str) -> bytes | None:
    """Read raw bytes from cache. Returns None if file doesn't exist."""
    path = _cache_path(cache_dir, source, filename)
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        return f.read()
