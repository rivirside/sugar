"""Detect duplicate compounds based on stereocenters within the same group."""


def _modifications_key(compound: dict) -> tuple:
    """Build a hashable key from a compound's modifications list."""
    mods = compound.get("modifications")
    if not mods:
        return ()
    return tuple(sorted((m["type"], m["position"]) for m in mods))


def check_duplicates(compounds: list[dict]) -> list[dict]:
    """Check for duplicate compounds in the compound list.

    Duplicates are identified as compounds within the same (type, carbons, modifications)
    group that have identical stereocenters.

    Returns a list of dicts describing each duplicate found.
    An empty list means no duplicates exist.
    """
    # Group by (type, carbons, modifications_key)
    groups: dict[tuple, list[dict]] = {}
    for c in compounds:
        ctype = c.get("type")
        carbons = c.get("carbons")
        mod_key = _modifications_key(c)
        key = (ctype, carbons, mod_key)
        groups.setdefault(key, []).append(c)

    duplicates = []
    for (ctype, carbons, _mod_key), group in groups.items():
        # Track seen stereocenters within the group
        seen: dict[tuple, dict] = {}
        for c in group:
            stereocenters_key = tuple(c.get("stereocenters", []))
            if stereocenters_key in seen:
                duplicates.append({
                    "original": seen[stereocenters_key],
                    "duplicate": c,
                    "type": ctype,
                    "carbons": carbons,
                    "stereocenters": list(stereocenters_key),
                })
            else:
                seen[stereocenters_key] = c

    return duplicates
