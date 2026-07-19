"""Sample-sheet driven input and genome-QC metadata import.

A sample sheet is a TSV with at least ``sample_id`` and ``annotation_file``,
optionally ``input_format`` and any number of metadata columns (group,
completeness, contamination, taxonomy, ...). It is the reproducible alternative
to scanning a directory.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .userdata import load_user_kos


def read_sample_sheet(path: str | Path) -> pd.DataFrame:
    """Read and validate a sample sheet (TSV). Returns a DataFrame."""
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "sample_id" not in df.columns:
        raise ValueError("sample sheet must have a 'sample_id' column")
    df["sample_id"] = df["sample_id"].str.strip()
    blank = df["sample_id"] == ""
    if blank.any():
        rows = [i + 2 for i in df.index[blank]]  # 1-indexed + header line
        raise ValueError(f"sample sheet has an empty sample_id at row(s): {rows}")
    if df["sample_id"].duplicated().any():
        dups = df.loc[df["sample_id"].duplicated(), "sample_id"].tolist()
        raise ValueError(f"duplicate sample_id(s): {dups}")
    return df


def load_kos_from_sheet(
    sheet: pd.DataFrame, base_dir: str | Path = ".",
) -> dict[str, set[str]]:
    """Load ``{sample_id: KO set}`` from a sample sheet's annotation files."""
    if "annotation_file" not in sheet.columns:
        raise ValueError("sample sheet must have an 'annotation_file' column")
    base = Path(base_dir)
    org_kos: dict[str, set[str]] = {}
    for _, row in sheet.iterrows():
        sid = row["sample_id"]
        annotation_file = str(row["annotation_file"]).strip()
        if not annotation_file:
            # Path(base) / "" resolves to base itself -- silently loading
            # every file in base_dir as this one sample's annotations.
            raise ValueError(
                f"sample sheet row for sample_id={sid!r} has an empty "
                "annotation_file")
        fmt = row.get("input_format") or "auto"
        loaded = load_user_kos(base / annotation_file, fmt)
        # A per-sample annotation file is one sample; take the union of its KOs.
        kos: set[str] = set().union(*loaded.values()) if loaded else set()
        org_kos[sid] = kos
    return org_kos


def sheet_metadata(sheet: pd.DataFrame) -> pd.DataFrame:
    """Return the metadata columns of a sample sheet, indexed by sample_id."""
    meta = sheet.set_index("sample_id")
    drop = [c for c in ("annotation_file", "input_format") if c in meta.columns]
    return meta.drop(columns=drop)


def import_checkm2(path: str | Path) -> pd.DataFrame:
    """Parse a CheckM2 quality_report.tsv -> completeness / contamination."""
    df = pd.read_csv(path, sep="\t")
    name = _first_col(df, ["Name", "name", "genome", "Bin Id"])
    comp = _first_col(df, ["Completeness", "completeness"])
    cont = _first_col(df, ["Contamination", "contamination"])
    return pd.DataFrame({
        "completeness": pd.to_numeric(df[comp], errors="coerce").to_numpy(),
        "contamination": pd.to_numeric(df[cont], errors="coerce").to_numpy(),
    }, index=df[name].astype(str)).rename_axis("sample_id")


def import_gtdbtk(path: str | Path) -> pd.DataFrame:
    """Parse a GTDB-Tk summary.tsv -> a 'taxonomy' column keyed by genome."""
    df = pd.read_csv(path, sep="\t")
    name = _first_col(df, ["user_genome", "genome", "Name"])
    clf = _first_col(df, ["classification", "Classification"])
    return pd.DataFrame({"taxonomy": df[clf].astype(str).to_numpy()},
                        index=df[name].astype(str)).rename_axis("sample_id")


def _first_col(df: pd.DataFrame, candidates: list[str]) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"expected one of {candidates} in columns {list(df.columns)}")
