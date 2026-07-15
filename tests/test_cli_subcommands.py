"""Offline tests for the post-processing subcommands (no network)."""
import json

import pandas as pd

import mpph.cli as cli


def _make_results(tmp_path, slug="Demo"):
    pd.DataFrame(
        {"00010": [1, 1, 1, 0], "00020": [1, 0, 0, 0], "00030": [1, 1, 0, 0]},
        index=["oA", "oB", "oC", "oD"],
    ).to_csv(tmp_path / f"{slug}_matrix.csv", index_label="organism")
    pd.DataFrame({"feature_id": ["00010", "00020", "00030"],
                  "name": ["a", "b", "c"],
                  "category": ["Carbohydrate metabolism"] * 3}).to_csv(
        tmp_path / f"{slug}_features.csv", index=False)
    pd.DataFrame({"organism": ["oA", "oB", "oC", "oD"],
                  "n_annotated_features": [3, 2, 1, 0],
                  "status": ["included"] * 4}).to_csv(
        tmp_path / f"{slug}_qc.csv", index=False)
    (tmp_path / f"{slug}_manifest.json").write_text(
        json.dumps({"mode": "presence", "mpph_version": "test"}))
    return slug


def test_pan_subcommand(tmp_path):
    _make_results(tmp_path)
    assert cli.main(["pan", "--results", str(tmp_path)]) == 0
    classes = pd.read_csv(tmp_path / "Demo_pan_classes.csv")
    assert set(classes["pan_class"]) <= {"core", "soft-core", "shell", "cloud"}
    assert (tmp_path / "Demo_accumulation.csv").exists()


def test_report_subcommand(tmp_path):
    _make_results(tmp_path)
    assert cli.main(["report", "--results", str(tmp_path)]) == 0
    assert (tmp_path / "Demo_report.html").exists()


def test_compare_subcommand(tmp_path):
    _make_results(tmp_path)
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample_id\tgroup\noA\tX\noB\tX\noC\tY\noD\tY\n")
    rc = cli.main(["compare", "--results", str(tmp_path), "--metadata", str(meta),
                   "--group-column", "group", "--group-a", "X", "--group-b", "Y"])
    assert rc == 0
    out = list(tmp_path.glob("Demo_differential_*.csv"))
    assert out and "q_value" in pd.read_csv(out[0]).columns


def test_validate_subcommand(tmp_path):
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample_id\tannotation_file\nA\ta.txt\nB\tb.txt\n")
    assert cli.main(["validate", "--samples", str(sheet)]) == 0


def test_backward_compatible_shortcut_routes_to_run(tmp_path, monkeypatch):
    # `mpph <taxon>` must still reach run(); stub run to avoid network.
    called = {}

    def _stub(args):
        called["taxon"] = args.taxon
        return 0

    monkeypatch.setattr(cli, "run", _stub)
    assert cli.main(["Prochlorococcus", "--completeness"]) == 0
    assert called["taxon"] == "Prochlorococcus"
