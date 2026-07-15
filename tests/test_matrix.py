"""Tests for matrix construction and filtering."""
from mpph.matrix import (
    build_completeness_matrix,
    build_presence_matrix,
    filter_matrix,
)


def test_presence_matrix_shape_and_values():
    org_pathways = {
        "org A": {"00010": "Glycolysis", "00020": "TCA"},
        "org B": {"00010": "Glycolysis"},
    }
    df, names = build_presence_matrix(org_pathways)
    assert df.shape == (2, 2)
    assert df.loc["org A", "00020"] == 1
    assert df.loc["org B", "00020"] == 0
    assert names["00010"] == "Glycolysis"


def test_completeness_matrix_values():
    org_kos = {"m1": {"K00001", "K00002"}, "m2": {"K00001"}}
    module_defs = {"M00001": ("K00001 K00002", "Carbohydrate metabolism")}
    df, _ = build_completeness_matrix(org_kos, module_defs)
    assert df.loc["m1", "M00001"] == 1.0
    assert df.loc["m2", "M00001"] == 0.5


def test_filter_drop_core():
    org_pathways = {
        "a": {"00010": "x", "00020": "y"},
        "b": {"00010": "x"},
    }
    df, _ = build_presence_matrix(org_pathways)
    filtered = filter_matrix(df, drop_core=True)
    assert "00010" not in filtered.columns  # present in all -> dropped
    assert "00020" in filtered.columns


def test_filter_min_prevalence():
    org_pathways = {
        "a": {"00010": "x", "00020": "y"},
        "b": {"00010": "x"},
        "c": {},
    }
    df, _ = build_presence_matrix(org_pathways)
    filtered = filter_matrix(df, min_prevalence=0.5)
    assert "00010" in filtered.columns   # 2/3
    assert "00020" not in filtered.columns  # 1/3


def test_filter_prevalence_state_complete():
    import pandas as pd
    # m1 fully complete in all; m2 only partially complete in most
    df = pd.DataFrame({"m1": [1.0, 1.0, 1.0], "m2": [0.5, 0.5, 1.0]},
                      index=list("abc"))
    # "any": m2 is detectable (>0) in all -> prevalence 1.0, kept
    assert "m2" in filter_matrix(df, min_prevalence=0.5,
                                 prevalence_state="any").columns
    # "complete": m2 fully complete in only 1/3 -> dropped; m1 kept
    kept = filter_matrix(df, min_prevalence=0.5, prevalence_state="complete")
    assert "m1" in kept.columns
    assert "m2" not in kept.columns
