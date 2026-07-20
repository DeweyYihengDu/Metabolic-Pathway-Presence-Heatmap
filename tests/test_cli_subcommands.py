"""Offline tests for the post-processing subcommands (no network)."""
import json
from pathlib import Path

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


def test_report_subcommand_max_cells_flag_reaches_build_report(tmp_path, monkeypatch):
    _make_results(tmp_path)
    seen = {}

    def _stub(outdir, slug, *, max_cells):
        seen["max_cells"] = max_cells
        return Path(outdir) / f"{slug}_report.html"

    import mpph.report as report_mod
    monkeypatch.setattr(report_mod, "build_report", _stub)
    cli.main(["report", "--results", str(tmp_path), "--max-cells", "123"])
    assert seen["max_cells"] == 123


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
    (tmp_path / "a.txt").write_text("K00001\n")
    (tmp_path / "b.txt").write_text("K00002\n")
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample_id\tannotation_file\nA\ta.txt\nB\tb.txt\n")
    assert cli.main(["validate", "--samples", str(sheet)]) == 0


def test_validate_subcommand_reports_missing_annotation_file(tmp_path, capsys):
    (tmp_path / "a.txt").write_text("K00001\n")
    sheet = tmp_path / "s.tsv"
    # B's annotation_file is never created.
    sheet.write_text("sample_id\tannotation_file\nA\ta.txt\nB\tmissing.txt\n")
    rc = cli.main(["validate", "--samples", str(sheet)])
    assert rc == 1
    assert "missing.txt" in capsys.readouterr().err


def test_validate_subcommand_reports_empty_annotation_file(tmp_path, capsys):
    (tmp_path / "a.txt").write_text("K00001\n")
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample_id\tannotation_file\nA\ta.txt\nB\t\n")
    rc = cli.main(["validate", "--samples", str(sheet)])
    assert rc == 1
    assert "empty annotation_file" in capsys.readouterr().err


def test_apply_genome_completeness_qc_drops_below_threshold():
    df = pd.DataFrame({"f1": [1.0, 1.0, 1.0]}, index=["mag1", "mag2", "mag3"])
    qc_meta = pd.DataFrame({"completeness": [98.5, 45.0, 85.0]},
                           index=["mag1", "mag2", "mag3"])
    out, excluded, low_included = cli.apply_genome_completeness_qc(
        df, qc_meta, min_completeness=50.0)
    assert list(out.index) == ["mag1", "mag3"]
    assert excluded == ["mag2"]
    assert low_included == ["mag3"]  # included but still <90%


def test_apply_genome_completeness_qc_no_filter_still_warns():
    df = pd.DataFrame({"f1": [1.0, 1.0]}, index=["mag1", "mag3"])
    qc_meta = pd.DataFrame({"completeness": [98.5, 85.0]}, index=["mag1", "mag3"])
    out, excluded, low_included = cli.apply_genome_completeness_qc(
        df, qc_meta, min_completeness=None)
    assert list(out.index) == ["mag1", "mag3"]  # nothing dropped
    assert excluded == []
    assert low_included == ["mag3"]


def test_apply_genome_completeness_qc_without_metadata_is_a_no_op():
    df = pd.DataFrame({"f1": [1.0]}, index=["mag1"])
    out, excluded, low_included = cli.apply_genome_completeness_qc(
        df, None, min_completeness=50.0)
    assert out is df
    assert excluded == [] and low_included == []


def test_run_argparser_accepts_qc_flags():
    args = cli.build_parser().parse_args([
        "run", "Vibrio", "--completeness", "--qc-metadata", "qc.tsv",
        "--min-genome-completeness", "70",
    ])
    assert args.qc_metadata == "qc.tsv"
    assert args.min_genome_completeness == 70.0


def test_run_argparser_qc_flags_default_to_none():
    args = cli.build_parser().parse_args(["run", "Vibrio"])
    assert args.qc_metadata is None
    assert args.min_genome_completeness is None


def test_traits_subcommand_writes_manifest_with_panel_provenance(tmp_path):
    users = tmp_path / "genomes"
    users.mkdir()
    (users / "mag1.txt").write_text("K02586\tK02588\tK02591\n")  # nitrogen_fixation
    outdir = tmp_path / "out"
    rc = cli.main(["traits", "--user", str(users), "--panel", "biogeochemistry",
                   "--outdir", str(outdir)])
    assert rc == 0
    manifest = json.loads((outdir / "biogeochemistry_manifest.json").read_text())
    assert manifest["kegg_release"] == "n/a (user data)"
    assert manifest["n_organisms"] == 1
    assert manifest["panel"]["panel_version"] == "1.0.0"
    assert len(manifest["panel"]["panel_sha256"]) == 64  # sha256 hex digest
    assert manifest["panel"]["panel_path"].endswith("biogeochemistry.json")
    assert "numpy" in manifest["dependency_versions"]
    assert "git_commit" in manifest  # present (value is None outside a checkout)


def test_run_provenance_reports_dependency_versions():
    prov = cli._run_provenance()
    deps = prov["dependency_versions"]
    assert set(deps) == {"numpy", "pandas", "scipy", "requests"}
    assert all(isinstance(v, str) and v for v in deps.values())
    assert "git_commit" in prov  # None or a str; never a missing key


def test_run_provenance_git_commit_is_none_outside_a_checkout(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)  # an empty tmp_path is never a git checkout
    assert cli._run_provenance()["git_commit"] is None


def test_backward_compatible_shortcut_routes_to_run(tmp_path, monkeypatch):
    # `mpph <taxon>` must still reach run(); stub run to avoid network.
    called = {}

    def _stub(args):
        called["taxon"] = args.taxon
        return 0

    monkeypatch.setattr(cli, "run", _stub)
    assert cli.main(["Prochlorococcus", "--completeness"]) == 0
    assert called["taxon"] == "Prochlorococcus"
