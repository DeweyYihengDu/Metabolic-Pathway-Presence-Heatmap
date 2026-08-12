"""Verifies compare_annotation.py's precision/recall/agreement logic against
hand-constructed synthetic tool outputs with known expected values, before
trusting it on real KofamScan/eggNOG-mapper/mpph annotate output. Not part of
the installed mpph package -- run explicitly via `pytest benchmarks/`, not
picked up by the main `pytest -q` (see pyproject.toml's testpaths).
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from compare_annotation import (
    _strip_version,
    build_report,
    load_emapper_annotations,
    load_ground_truth,
    load_mapper_tsv,
    load_protein_id_map,
    pairwise_agreement,
    precision_recall_f1,
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_load_mapper_tsv(tmp_path):
    p = _write(tmp_path / "m.tsv", "gene1\tK00001\ngene2\tK00002\ngene2\tK00003\n")
    assert load_mapper_tsv(p) == {"gene1": {"K00001"}, "gene2": {"K00002", "K00003"}}


def test_load_mapper_tsv_skips_kofamscan_no_hit_lines(tmp_path):
    # Real KofamScan --format mapper output (confirmed against actual output
    # on eco/bsu/mja): one line per query gene even with zero significant
    # hits -- just the bare gene id, no tab, no KO. ~23% of a real genome's
    # lines look like this; must be skipped, not treated as a parse error.
    p = _write(tmp_path / "kofamscan.tsv",
              "gene1\tK00001\ngene_no_hit\ngene2\tK00002\n")
    result = load_mapper_tsv(p)
    assert result == {"gene1": {"K00001"}, "gene2": {"K00002"}}
    assert "gene_no_hit" not in result


def test_load_emapper_annotations_handles_multi_ko_and_dash(tmp_path):
    p = _write(tmp_path / "e.annotations", (
        "## eggNOG-mapper\n"
        "#query\tseed_ortholog\tevalue\tscore\tKEGG_ko\n"
        "gene1\tx\t1e-50\t100\tko:K00001\n"
        "gene2\tx\t1e-40\t90\tko:K00002,ko:K00003\n"
        "gene4\tx\t1e-30\t80\t-\n"
    ))
    result = load_emapper_annotations(p)
    assert result == {"gene1": {"K00001"}, "gene2": {"K00002", "K00003"}}
    assert "gene4" not in result  # "-" means no call, not an empty-set call


def test_load_ground_truth_strips_org_and_ko_prefixes(tmp_path):
    p = _write(tmp_path / "t.tsv", "org:gene1\tko:K00001\norg:gene2\tko:K00002\n")
    assert load_ground_truth(p, "org") == {"gene1": {"K00001"}, "gene2": {"K00002"}}


def test_strip_version():
    assert _strip_version("NP_414542.1") == "NP_414542"
    assert _strip_version("NP_414542.12") == "NP_414542"
    assert _strip_version("NP_414542") == "NP_414542"  # no version suffix -- unchanged


def test_load_protein_id_map(tmp_path):
    # Real shape of KEGG's conv/<org>/ncbi-proteinid.
    p = _write(tmp_path / "map.tsv",
              "ncbi-proteinid:NP_414542\teco:b0001\nncbi-proteinid:NP_414543\teco:b0002\n")
    assert load_protein_id_map(p, "eco") == {"b0001": "NP_414542", "b0002": "NP_414543"}


def test_load_ground_truth_translates_via_protein_id_map(tmp_path):
    # This is the exact real-world bug this test guards against: KEGG's
    # link/ko/<org> ground truth is keyed by KEGG gene ids (b0001), but every
    # annotation tool here is keyed by whatever the genome FASTA's headers
    # were -- NCBI RefSeq protein accessions (NP_414542.1) -- a different id
    # namespace for the same gene, not a real accuracy difference. Without
    # --protein-id-map, comparing the two directly silently produces zero
    # true positives for every tool regardless of actual correctness.
    truth_path = _write(tmp_path / "truth.tsv", "eco:b0001\tko:K00001\neco:b0002\tko:K00002\n")
    id_map = {"b0001": "NP_414542", "b0002": "NP_414543"}
    result = load_ground_truth(truth_path, "eco", protein_id_map=id_map)
    assert result == {"NP_414542": {"K00001"}, "NP_414543": {"K00002"}}


def test_load_ground_truth_drops_genes_with_no_protein_id(tmp_path):
    # e.g. RNA genes with no protein product -- can't be matched against a
    # protein FASTA at all, so dropping (not keeping under an untranslatable
    # KEGG id that would never match anything) is the correct behaviour.
    truth_path = _write(tmp_path / "truth.tsv",
                        "eco:b0001\tko:K00001\neco:rna_gene\tko:K00002\n")
    id_map = {"b0001": "NP_414542"}  # rna_gene deliberately absent
    result = load_ground_truth(truth_path, "eco", protein_id_map=id_map)
    assert result == {"NP_414542": {"K00001"}}


def test_precision_recall_f1_known_values():
    # 2 correct, 0 wrong, 2 missed -> precision 1.0, recall 0.5, F1 2/3.
    predicted = {"gene1": {"K00001"}, "gene2": {"K00002"}}
    truth = {"gene1": {"K00001"}, "gene2": {"K00002", "K00003"}, "gene4": {"K00004"}}
    result = precision_recall_f1(predicted, truth)
    assert result["tp"] == 2 and result["fp"] == 0 and result["fn"] == 2
    assert result["precision"] == 1.0
    assert result["recall"] == 0.5
    assert result["f1"] == pytest.approx(2 / 3)


def test_precision_recall_f1_with_a_false_positive():
    predicted = {"gene1": {"K00001"}, "gene3": {"K99999"}}  # K99999 is wrong
    truth = {"gene1": {"K00001"}}
    result = precision_recall_f1(predicted, truth)
    assert result["tp"] == 1 and result["fp"] == 1 and result["fn"] == 0
    assert result["precision"] == 0.5
    assert result["recall"] == 1.0


def test_pairwise_agreement_jaccard():
    a = {"gene1": {"K00001"}, "gene2": {"K00002"}}          # 2 pairs
    b = {"gene1": {"K00001"}, "gene2": {"K00002", "K00003"}, "gene3": {"K99999"}}  # 4 pairs
    result = pairwise_agreement(a, b)
    # shared = 2 (gene1/K00001, gene2/K00002); union = 4 (a is a subset of b)
    assert result == {"n_a": 2, "n_b": 4, "n_shared": 2, "jaccard": 0.5}


def test_pairwise_agreement_disjoint_is_zero():
    a = {"gene1": {"K00001"}}
    b = {"gene2": {"K00002"}}
    assert pairwise_agreement(a, b)["jaccard"] == 0.0


def test_build_report_with_ground_truth_includes_precision_recall():
    calls = {"gene1": {"K00001"}}
    truth = {"gene1": {"K00001"}}
    accuracy, agreement = build_report("org", 10, calls, calls, calls, truth)
    assert "precision" in accuracy.columns
    assert accuracy["precision"].iloc[0] == 1.0
    assert len(agreement) == 3  # 3 tool pairs


def test_build_report_without_ground_truth_omits_precision_recall():
    # A MAG has no KEGG-curated KO assignments to score against. Coverage and
    # tool agreement are still meaningful; precision/recall are not, and must
    # be absent rather than silently reported as 0 (which would read as "the
    # tool got everything wrong" instead of "there is nothing to compare to").
    calls = {"gene1": {"K00001"}, "gene2": {"K00002"}}
    accuracy, agreement = build_report("mag", 10, calls, calls, calls, None)
    assert "precision" not in accuracy.columns
    assert "recall" not in accuracy.columns
    assert "tp" not in accuracy.columns
    assert accuracy["coverage"].iloc[0] == 0.2  # 2 called / 10 genes
    assert accuracy["n_calls"].iloc[0] == 2
    assert len(agreement) == 3  # agreement still computed
    assert agreement["jaccard"].iloc[0] == 1.0  # identical inputs
