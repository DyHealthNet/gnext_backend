# backend/utils/pheno_cache.py
from __future__ import annotations
import os
import threading
import pandas as pd
from typing import Dict, Tuple

# Shared state per process
_LOCK = threading.Lock()
_CACHE: Dict[Tuple[str, float], pd.DataFrame] = {}
_LAST_KEY: Tuple[str, float] | None = None   # fast path for one PHENO_FILE

def get_pheno_df(path: str) -> pd.DataFrame:
    """Return a cached DataFrame; reloads if the file changed on disk."""
    global _LAST_KEY
    mtime = os.path.getmtime(path)
    key = (path, mtime)

    # Fast path
    if _LAST_KEY == key and key in _CACHE:
        return _CACHE[key]

    with _LOCK:
        # Check again inside the lock
        if key in _CACHE:
            _LAST_KEY = key
            return _CACHE[key]

        # Drop old versions of the same path (if any)
        dead_keys = [k for k in _CACHE if k[0] == path and k != key]
        for k in dead_keys:
            _CACHE.pop(k, None)

        df = pd.read_csv(path)
        _CACHE[key] = df
        _LAST_KEY = key
        return df

def get_pheno_map(path: str) -> dict[str, tuple[str, str]]:
    """
    Convenience: { phenocode -> (category, description) }
    """
    df = get_pheno_df(path)
    # Using .itertuples is ~2x faster than iterrows for mapping
    return {row.phenocode: (row.category, row.description) for row in df.itertuples(index=False)}