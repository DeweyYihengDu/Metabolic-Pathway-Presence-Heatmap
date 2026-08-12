"""Aggregate every per-organism comparison into the paper-ready tables.

Combines the individual `report/annotation_<label>.csv` accuracy tables and
`*_agreement.csv` files written by compare_annotation.py, and pulls wall-clock
time and peak memory out of the `/usr/bin/time -v` logs each tool run wrote.

Usage:
    python summarize.py --report-dir report --results-dir results --out-dir report
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

# Tier assignment -- kept here (not inferred) so the distinction that matters
# scientifically is explicit rather than guessed from a filename.
MODEL_ORGANISMS = {"eco", "bsu", "mja"}
NONMODEL_ORGANISMS = {"vbs", "lbac"}

_ELAPSED = re.compile(r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)")
_MAXRSS = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")


def _tier(label: str) -> str:
    if label in MODEL_ORGANISMS:
        return "1_model_organism"
    if label in NONMODEL_ORGANISMS:
        return "2_nonmodel_genome"
    return "3_mag"


def parse_time_log(path: Path) -> dict:
    """Wall-clock seconds and peak RSS (MB) from a `/usr/bin/time -v` log."""
    if not path.exists():
        return {}
    text = path.read_text(encoding="utf-8", errors="ignore")
    out: dict = {}
    if (m := _ELAPSED.search(text)):
        parts = [float(p) for p in m.group(1).split(":")]
        seconds = 0.0
        for part in parts:          # h:mm:ss or m:ss
            seconds = seconds * 60 + part
        out["wall_clock_s"] = round(seconds, 1)
    if (m := _MAXRSS.search(text)):
        out["peak_memory_mb"] = round(int(m.group(1)) / 1024, 1)
    return out


def collect_accuracy(report_dir: Path, results_dir: Path) -> pd.DataFrame:
    frames = []
    for csv in sorted(report_dir.glob("annotation_*.csv")):
        if csv.name.endswith("_agreement.csv"):
            continue
        df = pd.read_csv(csv)
        label = csv.stem.removeprefix("annotation_")
        df["tier"] = _tier(label)
        tool_to_prefix = {"mpph_annotate": "mpph", "kofamscan": "kofamscan",
                          "eggnog_mapper": "emapper"}
        timings = [parse_time_log(results_dir / f"{tool_to_prefix[t]}_{label}.time.log")
                   for t in df["tool"]]
        for key in ("wall_clock_s", "peak_memory_mb"):
            df[key] = [t.get(key) for t in timings]
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    lead = ["tier", "organism", "tool"]
    return out[lead + [c for c in out.columns if c not in lead]].sort_values(
        ["tier", "organism", "tool"]).reset_index(drop=True)


def collect_agreement(report_dir: Path) -> pd.DataFrame:
    frames = []
    for csv in sorted(report_dir.glob("annotation_*_agreement.csv")):
        df = pd.read_csv(csv)
        df["tier"] = _tier(csv.stem.removeprefix("annotation_").removesuffix("_agreement"))
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    lead = ["tier", "organism", "tool_a", "tool_b"]
    return out[lead + [c for c in out.columns if c not in lead]].sort_values(
        ["tier", "organism"]).reset_index(drop=True)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--report-dir", default="report")
    p.add_argument("--results-dir", default="results")
    p.add_argument("--out-dir", default="report")
    args = p.parse_args()

    report_dir, results_dir = Path(args.report_dir), Path(args.results_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    accuracy = collect_accuracy(report_dir, results_dir)
    agreement = collect_agreement(report_dir)
    accuracy.to_csv(out_dir / "summary_accuracy.csv", index=False)
    agreement.to_csv(out_dir / "summary_agreement.csv", index=False)

    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print("=== accuracy / coverage / runtime ===")
        print(accuracy.to_string(index=False))
        print("\n=== pairwise tool agreement ===")
        print(agreement.to_string(index=False))
    print(f"\nWrote {out_dir / 'summary_accuracy.csv'} and "
          f"{out_dir / 'summary_agreement.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
