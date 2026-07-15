"""Tests for figure labelling (rendered to SVG so text is inspectable)."""
import pandas as pd

from mpph.plot import (
    plot_accumulation,
    plot_matrix,
    plot_prevalence,
    plot_volcano,
)


def _matrix():
    return pd.DataFrame(
        {"M1": [1.0, 0.5, 0.0], "M2": [0.2, 1.0, 0.8], "M3": [0.0, 0.3, 1.0]},
        index=["orgA", "orgB", "orgC"],
    )


def test_custom_value_and_strip_labels(tmp_path):
    out = tmp_path / "fig.svg"
    plot_matrix(_matrix(), {"M1": "Nitrogen cycle", "M2": "Sulfur cycle",
                            "M3": "Carbon cycle"},
                out, "Traits", "sub", cluster=True, mode="completeness",
                value_label="trait completeness", strip_label="trait category")
    svg = out.read_text(encoding="utf-8")
    assert "trait completeness" in svg
    assert "trait category" in svg
    assert "module completeness" not in svg


def test_default_presence_label(tmp_path):
    out = tmp_path / "pres.svg"
    plot_matrix(pd.DataFrame({"00010": [1, 0], "00020": [1, 1]},
                             index=["a", "b"]),
                {"00010": "Carbohydrate metabolism",
                 "00020": "Energy metabolism"},
                out, "Presence", "sub", mode="presence")
    svg = out.read_text(encoding="utf-8")
    assert "KEGG functional category" in svg


def test_analysis_plots_produce_files(tmp_path):
    diff = pd.DataFrame({"feature_id": ["a", "b", "c"],
                         "prevalence_diff": [0.8, -0.2, 0.1],
                         "q_value": [0.001, 0.9, 0.4]})
    plot_volcano(diff, tmp_path / "v.png", "A vs B")
    assert (tmp_path / "v.png").exists()

    pan = pd.DataFrame({"feature_id": ["a", "b", "c"],
                        "prevalence": [1.0, 0.5, 0.1],
                        "pan_class": ["core", "shell", "cloud"]})
    plot_prevalence(pan, tmp_path / "p.png", "pan")
    assert (tmp_path / "p.png").exists()

    acc = pd.DataFrame({"n_genomes": [1, 2, 3], "pan_mean": [3, 5, 6],
                        "core_mean": [3, 2, 1]})
    plot_accumulation(acc, tmp_path / "a.png", "acc")
    assert (tmp_path / "a.png").exists()
