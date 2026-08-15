"""Verifies arbitrate.py's candidate rules against hand-constructed cases with
known expected winners, before trusting them on real KOfam output. Run via
`pytest benchmarks/`, not the main `pytest -q` (see pyproject's testpaths).

These test the *offline evaluation* of candidate rules. The rule that shipped
is tested separately against the package implementation in
`tests/test_kofam.py::test_arbitration_*`.
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from arbitrate import apply_rule, evaluate


def mk(rows):
    df = pd.DataFrame(rows, columns=["gene","ko","delta","threshold","is_true"])
    df["score"] = df["delta"] + df["threshold"]
    df["relative"] = df["score"]/df["threshold"]
    return df

def test_top_delta_keeps_one_per_gene():
    df = mk([("g1","K1",5.0,100.0,True),("g1","K2",1.0,100.0,False),
             ("g2","K3",3.0,100.0,True)])
    out = apply_rule(df,"top_delta",tie_window=5.0,min_gap=10.0)
    assert set(out.ko)=={"K1","K3"} and len(out)==2

def test_relative_and_absolute_disagree_when_thresholds_differ():
    # +10 bits on a threshold of 1000 is a weaker signal than +8 on 50.
    df = mk([("g1","K_big",10.0,1000.0,False),("g1","K_small",8.0,50.0,True)])
    assert apply_rule(df,"top_delta",tie_window=5.0,min_gap=10.0).ko.iloc[0]=="K_big"
    assert apply_rule(df,"top_relative",tie_window=5.0,min_gap=10.0).ko.iloc[0]=="K_small"

def test_tie_window_keeps_near_ties():
    df = mk([("g1","K1",10.0,100.0,True),("g1","K2",7.0,100.0,False),
             ("g1","K3",1.0,100.0,False)])
    out = apply_rule(df,"top_delta_tie",tie_window=5.0,min_gap=10.0)
    assert set(out.ko)=={"K1","K2"}

def test_gap_gated_drops_ambiguous_genes_entirely():
    df = mk([("g1","K1",10.0,100.0,True),("g1","K2",9.0,100.0,False),  # gap 1 -> drop
             ("g2","K3",30.0,100.0,True),("g2","K4",5.0,100.0,False)]) # gap 25 -> keep
    out = apply_rule(df,"gap_gated",tie_window=5.0,min_gap=10.0)
    assert list(out.ko)==["K3"]

def test_gap_gated_keeps_unambiguous_single_call_genes():
    # A gene with only one passing KO has no runner-up and must be kept.
    df = mk([("g1","K1",2.0,100.0,True)])
    out = apply_rule(df,"gap_gated",tie_window=5.0,min_gap=10.0)
    assert list(out.ko)==["K1"]

def test_evaluate_counts_missed_truth_as_fn():
    df = mk([("g1","K1",5.0,100.0,True)])
    r = evaluate(df, n_true_pairs=3)
    assert r["tp"]==1 and r["fp"]==0 and r["fn"]==2
    assert r["precision"]==1.0 and r["recall"]==pytest.approx(1/3)

def test_all_rule_is_the_identity():
    df = mk([("g1","K1",5.0,100.0,True),("g1","K2",1.0,100.0,False)])
    assert len(apply_rule(df,"all",tie_window=5.0,min_gap=10.0))==2


@pytest.mark.parametrize("rule", ["top_delta","top_relative","top_delta_tie","gap_gated"])
def test_rules_are_invariant_to_input_row_order(rule):
    # Regression: `best` was computed on a sorted copy and compared back
    # against the unsorted frame, which raises on real data and would have
    # silently mis-selected if it had not. Every earlier test happened to
    # pass already-sorted input, so none of them caught it.
    rows = [("g1","K1",10.0,100.0,True),("g1","K2",2.0,100.0,False),
            ("g2","K3",30.0,50.0,True),("g2","K4",4.0,900.0,False),
            ("g3","K5",7.0,100.0,True)]
    sorted_in = mk(rows)
    shuffled = mk([rows[i] for i in (3,0,4,2,1)])
    a = apply_rule(sorted_in, rule, tie_window=5.0, min_gap=10.0)
    b = apply_rule(shuffled, rule, tie_window=5.0, min_gap=10.0)
    assert sorted(a.ko) == sorted(b.ko)
