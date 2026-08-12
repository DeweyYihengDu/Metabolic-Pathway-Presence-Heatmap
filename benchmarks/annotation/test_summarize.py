"""Verifies summarize.py's time-log parsing and tier assignment against
hand-constructed inputs with known expected values. Not part of the
installed mpph package -- run via `pytest benchmarks/`.
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from summarize import _tier, collect_accuracy, parse_time_log


def test_tier_assignment():
    assert _tier("eco") == "1_model_organism"
    assert _tier("vbs") == "2_nonmodel_genome"
    # anything not explicitly listed is a MAG -- the MAG labels are long
    # binning-tool filenames, not a fixed enumerable set.
    assert _tier("SRR12479784_metabat2_bin17") == "3_mag"


# Expected values are post-rounding: parse_time_log reports to 0.1 s
# deliberately (sub-decisecond precision is meaningless for multi-minute
# runs), so 2:49.25 -> 169.2, not 169.25.
@pytest.mark.parametrize("elapsed,expected_s", [
    ("2:49.25", 169.2),       # m:ss.ss -- the real format for these runs
    ("1:02:03", 3723.0),      # h:mm:ss
    ("0:45.00", 45.0),
])
def test_parse_time_log_elapsed_formats(tmp_path, elapsed, expected_s):
    p = tmp_path / "t.log"
    p.write_text(f"\tElapsed (wall clock) time (h:mm:ss or m:ss): {elapsed}\n"
                 f"\tMaximum resident set size (kbytes): 3826188\n")
    result = parse_time_log(p)
    assert result["wall_clock_s"] == expected_s
    assert result["peak_memory_mb"] == pytest.approx(3736.5, abs=0.1)


def test_parse_time_log_missing_file_is_empty_not_an_error(tmp_path):
    assert parse_time_log(tmp_path / "does_not_exist.log") == {}


def test_collect_accuracy_joins_timings_and_tiers(tmp_path):
    report, results = tmp_path / "report", tmp_path / "results"
    report.mkdir()
    results.mkdir()
    pd.DataFrame({
        "organism": ["eco", "eco"], "tool": ["mpph_annotate", "kofamscan"],
        "coverage": [0.76, 0.76], "n_calls": [3516, 3517],
    }).to_csv(report / "annotation_eco.csv", index=False)
    (results / "mpph_eco.time.log").write_text(
        "\tElapsed (wall clock) time (h:mm:ss or m:ss): 2:49.25\n"
        "\tMaximum resident set size (kbytes): 1024000\n")
    # kofamscan's log deliberately absent -- must yield NaN, not raise.

    out = collect_accuracy(report, results)
    assert list(out["tier"]) == ["1_model_organism"] * 2
    mpph_row = out[out["tool"] == "mpph_annotate"].iloc[0]
    assert mpph_row["wall_clock_s"] == 169.2  # rounded to 0.1 s, see above
    assert mpph_row["peak_memory_mb"] == pytest.approx(1000.0)
    assert pd.isna(out[out["tool"] == "kofamscan"].iloc[0]["wall_clock_s"])


def test_collect_accuracy_skips_agreement_files(tmp_path):
    report, results = tmp_path / "report", tmp_path / "results"
    report.mkdir()
    results.mkdir()
    pd.DataFrame({"organism": ["eco"], "tool": ["mpph_annotate"],
                  "coverage": [0.5], "n_calls": [1]}).to_csv(
        report / "annotation_eco.csv", index=False)
    # Same glob prefix -- must not be mistaken for an accuracy table.
    pd.DataFrame({"organism": ["eco"], "tool_a": ["mpph_annotate"],
                  "tool_b": ["kofamscan"], "jaccard": [1.0]}).to_csv(
        report / "annotation_eco_agreement.csv", index=False)

    out = collect_accuracy(report, results)
    assert len(out) == 1
    assert "jaccard" not in out.columns
