"""Tests for figure labelling (rendered to SVG so text is inspectable)."""
import pandas as pd

from mpph.plot import plot_matrix


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
