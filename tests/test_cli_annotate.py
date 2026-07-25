"""Tests for `mpph annotate`'s CLI wiring -- argument parsing, the
--setup-db path, the friendly missing-pyhmmer message, and manifest writing.
These monkeypatch mpph.kofam's own functions rather than re-testing the
search/scoring logic itself (that's tests/test_kofam.py and
tests/test_kofam_pyhmmer.py's job); no real network call and no real
pyhmmer search happens here regardless of what's installed.
"""
import builtins
import json

import pytest

import mpph.cli as cli
import mpph.kofam as kofam


def test_annotate_argparser_fasta_and_setup_db_are_mutually_exclusive_and_required():
    parser = cli.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["annotate"])
    with pytest.raises(SystemExit):
        parser.parse_args(["annotate", "--fasta", "x.faa", "--setup-db", "d"])


def test_annotate_argparser_defaults():
    args = cli.build_parser().parse_args(["annotate", "--fasta", "x.faa"])
    assert args.kofam_db == cli.DEFAULT_KOFAM_DB_STR
    assert args.out is None
    assert args.cpus == 0
    assert args.ko_subset is None
    assert args.overwrite_db is False


def test_cmd_annotate_setup_db_calls_download_kofam_db(tmp_path, monkeypatch):
    calls = {}

    def _fake_download(dest_dir, *, overwrite=False, **kwargs):
        calls["dest_dir"] = dest_dir
        calls["overwrite"] = overwrite
        return dest_dir

    monkeypatch.setattr(kofam, "download_kofam_db", _fake_download)
    dest = tmp_path / "kofam_db"
    rc = cli.main(["annotate", "--setup-db", str(dest), "--overwrite-db"])
    assert rc == 0
    assert str(calls["dest_dir"]) == str(dest)
    assert calls["overwrite"] is True


def test_cmd_annotate_missing_pyhmmer_gives_friendly_message(tmp_path, monkeypatch, capsys):
    real_import = builtins.__import__

    def _fake_import(name, *args, **kwargs):
        if name == "pyhmmer":
            raise ImportError("simulated: pyhmmer not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _fake_import)

    rc = cli.main(["annotate", "--fasta", str(tmp_path / "does_not_matter.faa")])
    assert rc == 1
    assert "pip install mpph[annotate]" in capsys.readouterr().err


def test_cmd_annotate_missing_kofam_db_gives_clean_error_not_traceback(tmp_path, capsys):
    # cmd_annotate's own pyhmmer-availability probe runs before this code
    # path, so real pyhmmer must importable to reach it -- the `annotate`
    # extra is deliberately not part of `dev`/CI's install, so skip rather
    # than fail where it's absent (same reasoning as test_kofam_pyhmmer.py).
    pytest.importorskip("pyhmmer")
    fasta = tmp_path / "x.faa"
    fasta.write_text(">g\nMKT\n")
    rc = cli.main(["annotate", "--fasta", str(fasta),
                  "--kofam-db", str(tmp_path / "no_such_db")])
    assert rc == 1
    assert "annotate --setup-db" in capsys.readouterr().err


def test_cmd_annotate_happy_path_writes_output_and_manifest(tmp_path, monkeypatch):
    # annotate_fasta is monkeypatched below (no real search happens), but
    # cmd_annotate's own `import pyhmmer` availability probe still runs
    # first and must succeed to reach that code -- see importorskip note
    # on test_cmd_annotate_missing_kofam_db_gives_clean_error_not_traceback.
    pytest.importorskip("pyhmmer")
    fasta = tmp_path / "genome.faa"
    fasta.write_text(">gene1\nMKT\n>gene2\nMKT\n")
    (tmp_path / "kofam_db").mkdir()
    (tmp_path / "kofam_db" / "ko_list").write_text(
        "knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\t"
        "alen\tmlen\teff_nseq\tre/pos\tdefinition\n"
        "K00001\t100.0\tfull\tfull\t0.9\t10\t10\t50\t50\t5.0\t0.5\ttest ko\n")

    canned = [kofam.Assignment("gene1", "K00001", 150.0),
              kofam.Assignment("gene2", "K00001", 120.0)]
    monkeypatch.setattr(kofam, "annotate_fasta",
                        lambda *a, **k: canned)

    out = tmp_path / "genome_annotated.tsv"
    rc = cli.main(["annotate", "--fasta", str(fasta),
                  "--kofam-db", str(tmp_path / "kofam_db"), "--out", str(out)])
    assert rc == 0
    assert out.read_text(encoding="utf-8").splitlines() == [
        "gene1\tK00001", "gene2\tK00001"]

    manifest = json.loads(
        out.with_name("genome_annotated_manifest.json").read_text())
    assert manifest["n_assignments"] == 2
    assert manifest["n_genes_annotated"] == 2
    assert manifest["fasta"] == str(fasta)
    assert manifest["outputs"] == {"annotations": "genome_annotated.tsv"}
    assert "git_commit" in manifest  # provenance wired in like run/traits/pathmap


def test_cmd_annotate_default_out_path_is_fasta_stem_annotated(tmp_path, monkeypatch):
    # Same pyhmmer-availability-probe reasoning as the happy-path test above.
    pytest.importorskip("pyhmmer")
    # cmd_annotate resolves the default relative to the CWD, not the fasta's
    # own directory -- chdir into a disposable tmp_path to check this safely
    # (the documented default is <fasta stem>_annotated.tsv).
    monkeypatch.chdir(tmp_path)
    fasta = tmp_path / "myMAG.faa"
    fasta.write_text(">gene1\nMKT\n")
    (tmp_path / "kofam_db").mkdir()
    (tmp_path / "kofam_db" / "ko_list").write_text(
        "knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\t"
        "alen\tmlen\teff_nseq\tre/pos\tdefinition\n")
    monkeypatch.setattr(kofam, "annotate_fasta", lambda *a, **k: [])

    rc = cli.main(["annotate", "--fasta", str(fasta),
                  "--kofam-db", str(tmp_path / "kofam_db")])
    assert rc == 0
    assert (tmp_path / "myMAG_annotated.tsv").exists()
