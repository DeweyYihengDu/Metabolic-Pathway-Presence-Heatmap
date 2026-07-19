"""Load user-supplied KO annotations so MPPH can plot your own genomes/MAGs.

Accepted inputs (auto-detected):

* **A directory** -- one file per genome/MAG. Each file is scanned for KO ids
  (``K\\d{5}``), so KofamScan ``--format mapper`` output (``gene<TAB>KO``) or a
  bare list of KO ids both work. The sample name is the file stem.
* **A single 2-column table** -- ``sample<TAB>KO`` (long form); KO ids are
  extracted per row and grouped by the first column.
* **An eggNOG-mapper ``.annotations`` file** -- KO ids are read from the
  ``KEGG_ko`` column (values like ``ko:K00001``).

Returns ``{sample_name: {KO, ...}}``.
"""
from __future__ import annotations

import re
from pathlib import Path

_KO = re.compile(r"K\d{5}")


def _kos_from_text(text: str) -> set[str]:
    return set(_KO.findall(text))


def _kos_from_eggnog(text: str) -> set[str]:
    """Extract KOs from an eggNOG-mapper ``.annotations`` file's KEGG_ko column."""
    header: list[str] | None = None
    idx = None
    kos: set[str] = set()
    for line in text.splitlines():
        if line.startswith("##") or not line.strip():
            continue
        fields = line.lstrip("#").split("\t")
        if header is None:
            if "KEGG_ko" in fields:
                header = fields
                idx = header.index("KEGG_ko")
            continue
        if idx is not None and idx < len(fields):
            kos.update(_KO.findall(fields[idx]))
    return kos


def _load_long_table(path: Path) -> tuple[dict[str, set[str]], int]:
    """Parse a ``sample<TAB>...KO...`` table (one record per line).

    Returns ``(samples, n_data_rows)`` -- the row count excludes blank,
    comment and KO-less lines, so it can be compared against ``len(samples)``
    to tell a genuine long table (sample ids repeat across rows) from a file
    that merely happens to have >=2 tab-separated fields per line.
    """
    samples: dict[str, set[str]] = {}
    n_rows = 0
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        kos = _KO.findall(line)
        if not kos:
            continue
        n_rows += 1
        sample = fields[0].strip() or path.stem
        samples.setdefault(sample, set()).update(kos)
    return samples, n_rows


def load_user_kos(path: str | Path, fmt: str = "auto") -> dict[str, set[str]]:
    """Load KO sets keyed by sample name from ``path`` (dir or file).

    ``fmt`` is ``auto`` (default), ``ko-list`` (scan any ``K#####``, whole
    file/each directory entry is one sample), ``long`` (force ``sample<TAB>KO``
    parsing), or ``eggnog`` (read the ``KEGG_ko`` column of eggNOG-mapper
    output). For a directory, each file is one sample (named by its stem).

    A single two-column file is ambiguous by construction: it could be a long
    table (``sample<TAB>KO``, the same sample id repeating over many rows) or
    a one-sample-per-file annotation table (e.g. KofamScan ``--format mapper``,
    ``gene<TAB>KO``, every row a *different* gene). ``auto`` only treats it as
    a long table when at least one sample id genuinely repeats across rows;
    otherwise the whole file is one sample (named by its stem), matching this
    tool's primary one-file-per-genome use. Pass ``fmt="long"`` to force the
    long-table reading for a genuine multi-sample table where every sample
    happens to contribute exactly one row.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"user KO input not found: {p}")
    extract = _kos_from_eggnog if fmt == "eggnog" else _kos_from_text

    if p.is_dir():
        samples: dict[str, set[str]] = {}
        files = sorted(f for f in p.iterdir()
                       if f.is_file() and not f.name.startswith("."))
        if not files:
            raise ValueError(f"no files found in directory: {p}")
        for f in files:
            kos = extract(f.read_text(encoding="utf-8", errors="ignore"))
            if kos:
                if f.stem in samples:
                    raise ValueError(
                        f"duplicate sample name {f.stem!r} in {p}: multiple "
                        "files share this stem (e.g. .txt and .tsv), so the "
                        "sample name is ambiguous. Rename one of them.")
                samples[f.stem] = kos
        if not samples:
            raise ValueError(f"no KO ids (K#####) found in any file under {p}")
        return samples

    text = p.read_text(encoding="utf-8", errors="ignore")
    if fmt == "eggnog":
        kos = _kos_from_eggnog(text)
        if not kos:
            raise ValueError(f"no KEGG_ko entries found in {p}")
        return {p.stem: kos}
    if fmt == "long":
        table, _n_rows = _load_long_table(p)
        if not table:
            raise ValueError(f"no KO ids (K#####) found in {p}")
        return table
    if fmt == "auto":
        # long table (sample<TAB>KO) vs. whole-file-is-one-sample: only
        # trust the long-table reading when a sample id actually repeats
        # across rows -- otherwise a per-gene KofamScan file (one row per
        # gene, every gene id unique) would silently turn each gene into its
        # own "sample". See fmt="long" to force this reading regardless.
        first = next((ln for ln in text.splitlines() if ln.strip()), "")
        if len(first.split("\t")) >= 2 and _KO.search(text):
            table, n_rows = _load_long_table(p)
            if len(table) < n_rows:
                return table
    kos = _kos_from_text(text)
    if not kos:
        raise ValueError(f"no KO ids (K#####) found in {p}")
    return {p.stem: kos}
