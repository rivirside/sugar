"""BRENDA enzyme database importer. Uses SOAP API per EC number."""

import hashlib
import logging
import os

from dotenv import load_dotenv
from pipeline.import_.cache import read_cache, write_cache, is_cache_fresh

logger = logging.getLogger(__name__)
BRENDA_WSDL = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"


def load_brenda_credentials() -> tuple[str, str] | None:
    load_dotenv()
    email = os.getenv("BRENDA_EMAIL")
    password = os.getenv("BRENDA_PASSWORD")
    if not email or not password:
        logger.warning("BRENDA credentials not found in .env file")
        return None
    return email, password


def fetch_brenda_kinetics(ec_numbers: list[str], cache_dir: str, refresh: bool = False) -> dict:
    credentials = load_brenda_credentials()
    if not credentials:
        logger.warning("Skipping BRENDA import (no credentials)")
        return {}

    email, password = credentials
    password_hash = hashlib.sha256(password.encode()).hexdigest()

    results = {}
    for ec in ec_numbers:
        cache_file = f"{ec.replace('.', '_')}.json"
        if not refresh and is_cache_fresh(cache_dir, "brenda", cache_file):
            cached = read_cache(cache_dir, "brenda", cache_file)
            if cached:
                results[ec] = cached
                continue
        try:
            kinetics = _fetch_ec_kinetics(email, password_hash, ec)
            write_cache(cache_dir, "brenda", cache_file, kinetics)
            results[ec] = kinetics
        except Exception as e:
            logger.warning("BRENDA fetch failed for EC %s: %s", ec, e)
    return results


def _fetch_ec_kinetics(email: str, password_hash: str, ec_number: str) -> dict:
    try:
        from zeep import Client
        client = Client(BRENDA_WSDL)
        km_raw = []
        kcat_raw = []
        try:
            km_result = client.service.getKmValue(email, password_hash, f"ecNumber*{ec_number}")
            if km_result:
                km_raw = parse_brenda_km_data(km_result if isinstance(km_result, list) else [])
        except Exception as e:
            logger.debug("BRENDA Km fetch failed for %s: %s", ec_number, e)
        try:
            kcat_result = client.service.getTurnoverNumber(email, password_hash, f"ecNumber*{ec_number}")
            if kcat_result:
                kcat_raw = parse_brenda_kcat_data(kcat_result if isinstance(kcat_result, list) else [])
        except Exception as e:
            logger.debug("BRENDA kcat fetch failed for %s: %s", ec_number, e)
        return {"ec_number": ec_number, "km_entries": km_raw, "kcat_entries": kcat_raw}
    except ImportError:
        logger.warning("zeep not installed, skipping BRENDA SOAP")
        return {"ec_number": ec_number, "km_entries": [], "kcat_entries": []}


def parse_brenda_km_data(km_entries: list) -> list[dict]:
    results = []
    for entry in km_entries:
        if isinstance(entry, dict):
            results.append({"ec_number": entry.get("ecNumber", ""), "km_mm": entry.get("kmValue"), "substrate": entry.get("substrate", ""), "organism": entry.get("organism", "")})
    return results


def parse_brenda_kcat_data(kcat_entries: list) -> list[dict]:
    results = []
    for entry in kcat_entries:
        if isinstance(entry, dict):
            results.append({"ec_number": entry.get("ecNumber", ""), "kcat_sec": entry.get("turnoverNumber"), "substrate": entry.get("substrate", ""), "organism": entry.get("organism", "")})
    return results
