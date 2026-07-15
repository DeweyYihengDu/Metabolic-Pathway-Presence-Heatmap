"""Low-level KEGG REST access: pooled session, retries, disk cache."""
from __future__ import annotations

import re
import time
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter

try:  # urllib3 ships with requests; import path differs across versions
    from urllib3.util.retry import Retry
except ImportError:  # pragma: no cover
    from requests.packages.urllib3.util.retry import Retry  # type: ignore

KEGG_API_BASE = "https://rest.kegg.jp"
REQUEST_TIMEOUT = 30  # seconds
DEFAULT_CACHE = Path(".mpph_cache")


def make_session() -> requests.Session:
    """A requests session that retries transient failures with backoff."""
    retry = Retry(
        total=5,
        backoff_factor=0.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"}),
    )
    session = requests.Session()
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    from . import __version__  # local import avoids a circular import at load
    session.headers.update(
        {"User-Agent": f"MPPH/{__version__} (+github.com/DeweyYihengDu)"})
    return session


def _cache_path(cache_dir: Path, endpoint: str) -> Path:
    return cache_dir / (re.sub(r"[^0-9A-Za-z]+", "_", endpoint) + ".tsv")


def kegg_get(
    session: requests.Session,
    endpoint: str,
    cache_dir: Path | None,
    *,
    refresh: bool = False,
    max_age: float | None = None,
) -> str:
    """GET ``{KEGG_API_BASE}/{endpoint}`` with optional disk caching.

    ``refresh`` forces a re-fetch; ``max_age`` (seconds) treats a cache entry
    older than that as stale. Returns the raw response text.
    """
    cache_file = _cache_path(cache_dir, endpoint) if cache_dir is not None else None
    if cache_file is not None and cache_file.exists() and not refresh:
        fresh = max_age is None or (time.time() - cache_file.stat().st_mtime) < max_age
        if fresh:
            return cache_file.read_text(encoding="utf-8")

    resp = session.get(f"{KEGG_API_BASE}/{endpoint}", timeout=REQUEST_TIMEOUT)
    resp.raise_for_status()
    text = resp.text

    if cache_file is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(text, encoding="utf-8")
    return text


def kegg_release(
    session: requests.Session, cache_dir: Path | None, *, refresh: bool = False
) -> str:
    """Return the KEGG release string (from the ``info/kegg`` endpoint)."""
    try:
        text = kegg_get(session, "info/kegg", cache_dir, refresh=refresh)
    except requests.RequestException:
        return "unknown"
    for line in text.splitlines():
        if "Release" in line or "release" in line:
            return line.strip()
    # Fall back to the first dated line (e.g. "pathway  587  2026/07/14").
    for line in text.splitlines():
        if re.search(r"\d{4}/\d{2}/\d{2}", line):
            return line.strip()
    return "unknown"
