"""Verifies score_precision.py's binning, truth-joining and shuffle control
against hand-constructed data with known expected values, before trusting it
on real KOfam output. Run via `pytest benchmarks/`.
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from score_precision import annotate_truth, bin_precision, load_hits, shuffle_control


def _hits(rows):
    return pd.DataFrame(rows, columns=["gene", "ko", "delta"])


def test_annotate_truth_marks_only_pairs_present_in_the_reference():
    df = _hits([("g1", "K00001", 5.0), ("g1", "K99999", -3.0),
                ("g2", "K00002", -1.0)])
    truth = {"g1": {"K00001"}, "g2": {"K00002"}}
    out = annotate_truth(df, truth)
    assert list(out["is_true"]) == [True, False, True]


def test_annotate_truth_is_pairwise_not_per_gene():
    # A gene present in the reference does not make every KO called on it
    # correct -- the unit is the (gene, KO) pair. Getting this wrong would
    # inflate sub-threshold precision, which is exactly the number at stake.
    df = _hits([("g1", "K00001", -2.0), ("g1", "K00002", -2.0)])
    truth = {"g1": {"K00001"}}
    assert list(annotate_truth(df, truth)["is_true"]) == [True, False]


def test_bin_precision_known_values():
    df = _hits([("g1", "K1", -2.0), ("g2", "K2", -3.0),
                ("g3", "K3", -7.0), ("g4", "K4", 1.0)])
    truth = {"g1": {"K1"}, "g4": {"K4"}}          # one true in [-5,0), one in [0,5)
    out = bin_precision(annotate_truth(df, truth), width=5.0, lo=-10.0, hi=5.0)
    by_lo = {row.delta_lo: row for row in out.itertuples()}
    assert by_lo[-5.0].n_hits == 2 and by_lo[-5.0].precision == 0.5
    assert by_lo[-10.0].n_hits == 1 and by_lo[-10.0].precision == 0.0
    assert by_lo[0.0].n_hits == 1 and by_lo[0.0].precision == 1.0


def test_bin_precision_skips_empty_bins_rather_than_reporting_zero():
    # An empty bin has undefined precision; emitting 0.0 would read as "no
    # true hits here" instead of "no data here".
    df = _hits([("g1", "K1", 1.0)])
    out = bin_precision(annotate_truth(df, {"g1": {"K1"}}),
                        width=5.0, lo=-20.0, hi=5.0)
    assert len(out) == 1
    assert out.iloc[0]["delta_lo"] == 0.0


def test_shuffle_control_preserves_bin_sizes():
    # The control must differ from the observed only in the KO labelling --
    # if it also changed how many hits land in each bin it would not be a
    # like-for-like background.
    df = _hits([("g1", "K1", -2.0), ("g2", "K2", -3.0), ("g3", "K3", -7.0)])
    truth = {"g1": {"K1"}, "g2": {"K2"}}
    scored = annotate_truth(df, truth)
    obs = bin_precision(scored, width=5.0, lo=-10.0, hi=0.0)
    ctl = shuffle_control(scored, truth, seed=0, width=5.0, lo=-10.0, hi=0.0)
    assert list(obs["n_hits"]) == list(ctl["n_hits"])


def test_shuffle_control_destroys_a_perfect_signal():
    # With many distinct KOs, permuting labels should almost surely break the
    # gene->KO correspondence, so the control precision must drop well below
    # the observed 1.0.
    rows = [(f"g{i}", f"K{i:05d}", -2.0) for i in range(60)]
    truth = {f"g{i}": {f"K{i:05d}"} for i in range(60)}
    scored = annotate_truth(_hits(rows), truth)
    assert bin_precision(scored, width=5.0, lo=-5.0, hi=0.0)["precision"].iloc[0] == 1.0
    ctl = shuffle_control(scored, truth, seed=0, width=5.0, lo=-5.0, hi=0.0)
    assert ctl["precision"].iloc[0] < 0.2


def test_load_hits_strips_refseq_version_suffix(tmp_path):
    # Must match the annotation benchmark's id handling or the join against
    # KEGG's conv table silently yields zero true positives.
    p = tmp_path / "h.tsv"
    p.write_text("gene\tko\tscore\tthreshold\tdelta\tscore_type\n"
                 "NP_414542.1\tK00001\t100.0\t90.0\t10.0\tfull\n", encoding="utf-8")
    assert load_hits(p)["gene"].iloc[0] == "NP_414542"


@pytest.mark.parametrize("delta,expected_lo", [(-0.001, -5.0), (0.0, 0.0), (4.999, 0.0)])
def test_bin_edges_put_zero_in_the_kept_band(delta, expected_lo):
    # delta == 0 means score exactly equals threshold, which mpph and
    # KofamScan both *accept* (>= not >), so it belongs in the kept band.
    df = _hits([("g1", "K1", delta)])
    out = bin_precision(annotate_truth(df, {}), width=5.0, lo=-10.0, hi=10.0)
    assert out.iloc[0]["delta_lo"] == expected_lo


# --- summarize_audit.py -------------------------------------------------

from summarize_audit import LOW_MARGIN, summarize


def _scored(rows):
    df = pd.DataFrame(rows, columns=["delta", "is_true"])
    return df


def test_summarize_splits_accepted_and_discarded_at_zero():
    df = _scored([(0.0, True), (-0.001, False), (10.0, True)])
    row, _ = summarize(df, "x")
    assert row["n_accepted"] == 2          # delta == 0 is accepted (>=, not >)
    assert row["n_discarded_within_60bits"] == 1


def test_fp_enrichment_is_one_when_errors_are_spread_evenly():
    # Half the calls low-margin, half high, same error rate in both -> the
    # low-margin band carries exactly its share, enrichment 1.0. Guards
    # against a definition that reports enrichment for a uniform genome.
    rows = ([(5.0, True)] * 5 + [(5.0, False)] * 5
            + [(500.0, True)] * 5 + [(500.0, False)] * 5)
    row, _ = summarize(_scored(rows), "x")
    assert row["frac_calls_low_margin"] == pytest.approx(0.5)
    assert row["fp_enrichment_low_margin"] == pytest.approx(1.0)


def test_fp_enrichment_exceeds_one_when_errors_concentrate_near_threshold():
    rows = ([(5.0, False)] * 8 + [(5.0, True)] * 2
            + [(500.0, True)] * 88 + [(500.0, False)] * 2)
    row, _ = summarize(_scored(rows), "x")
    assert row["fp_enrichment_low_margin"] > 3.0
    assert row["precision_low_margin"] == pytest.approx(0.2)


def test_summarize_curve_fractions_sum_to_one():
    rows = [(1.0, True), (7.0, False), (15.0, True), (30.0, True),
            (50.0, False), (80.0, True), (150.0, True), (900.0, True)]
    _, curve = summarize(_scored(rows), "x")
    assert curve["frac_of_calls"].sum() == pytest.approx(1.0)
    assert curve["n_calls"].sum() == 8


def test_low_margin_band_matches_the_documented_constant():
    # The headline statistic depends on this cut; if it moves, the numbers
    # quoted in the write-up silently stop matching the code.
    assert LOW_MARGIN == 20.0


# --- calibration_transfer.py --------------------------------------------

import numpy as np
from calibration_transfer import brier, ece, fit_logistic, predict


def test_ece_is_zero_for_a_perfectly_calibrated_predictor():
    # Half the items at confidence 0.5 and truly right half the time.
    p = np.full(1000, 0.5)
    y = np.array([1.0, 0.0] * 500)
    assert ece(p, y) == pytest.approx(0.0, abs=1e-12)


def test_ece_detects_overconfidence():
    # Claims 0.95, right only 50% of the time -> ECE ~ 0.45.
    p = np.full(1000, 0.95)
    y = np.array([1.0, 0.0] * 500)
    assert ece(p, y) == pytest.approx(0.45, abs=0.01)


def test_brier_rewards_sharpness_a_constant_cannot_fake():
    y = np.array([1.0] * 500 + [0.0] * 500)
    vague = np.full(1000, 0.5)                       # calibrated but useless
    sharp = np.array([0.9] * 500 + [0.1] * 500)      # calibrated and informative
    assert brier(sharp, y) < brier(vague, y)


def test_fit_logistic_recovers_a_monotone_increasing_relationship():
    rng = np.random.default_rng(0)
    delta = rng.uniform(0, 500, 4000)
    # Truth: probability rises with margin.
    p_true = 1.0 / (1.0 + np.exp(-(-2.0 + 0.9 * np.log1p(delta))))
    y = (rng.uniform(size=delta.size) < p_true).astype(float)
    a, b = fit_logistic(delta, y)
    assert b > 0                                    # increasing in margin
    fitted = predict(np.array([1.0, 500.0]), a, b)
    assert fitted[1] > fitted[0]
    assert abs(fitted[1] - p_true[np.argmax(delta)]) < 0.15


def test_predict_is_bounded_and_monotone():
    a, b = -2.0, 0.9
    p = predict(np.array([0.0, 10.0, 100.0, 10000.0]), a, b)
    assert np.all((p >= 0) & (p <= 1))
    assert np.all(np.diff(p) > 0)


def test_predict_clips_negative_margins_rather_than_producing_nan():
    # Only accepted calls are calibrated, but a negative margin must not
    # produce a NaN through log1p of a negative number.
    assert np.isfinite(predict(np.array([-5.0]), -2.0, 0.9)).all()


from calibration_transfer import bootstrap_ci


def test_bootstrap_ci_brackets_zero_when_the_two_predictors_are_identical():
    rng = np.random.default_rng(0)
    y = (rng.uniform(size=400) < 0.8).astype(float)
    p = np.full(400, 0.8)
    ci = bootstrap_ci(p, p.copy(), y, n_boot=300, seed=1)
    assert ci["ece_improvement_lo"] <= 0.0 <= ci["ece_improvement_hi"]
    assert ci["brier_improvement_lo"] == pytest.approx(0.0, abs=1e-12)


def test_bootstrap_ci_excludes_zero_for_a_clearly_better_predictor():
    y = np.array([1.0] * 500 + [0.0] * 500)
    good = np.array([0.95] * 500 + [0.05] * 500)   # informative
    const = np.full(1000, 0.5)                     # uninformative
    ci = bootstrap_ci(good, const, y, n_boot=300, seed=1)
    assert ci["brier_improvement_lo"] > 0.0


# --- neighbour_calibrator.py: hand-rolled PAVA -------------------------------

from neighbour_calibrator import (
    apply_platt,
    isotonic_on_neighbours,
    murphy_decomposition,
)


def test_isotonic_output_is_monotone_non_decreasing():
    # The defining property. A violation would silently make a "calibrated"
    # score rank calls worse than the raw margin it was built from.
    rng = np.random.default_rng(0)
    p = np.sort(rng.uniform(size=200))
    y = (rng.uniform(size=200) < p).astype(float)
    out = isotonic_on_neighbours(p, y, p)
    assert np.all(np.diff(out) >= -1e-12)


def test_isotonic_pools_adjacent_violators_to_their_mean():
    # Hand-checked: y = 0,1,0,1,1 at increasing x pools to 0, .5, .5, 1, 1.
    p = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    y = np.array([0.0, 1.0, 0.0, 1.0, 1.0])
    out = isotonic_on_neighbours(p, y, p)
    assert out == pytest.approx([0.0, 0.5, 0.5, 1.0, 1.0])


def test_isotonic_is_the_identity_on_already_calibrated_input():
    p = np.array([0.1, 0.3, 0.6, 0.9])
    y = np.array([0.0, 0.0, 1.0, 1.0])
    assert isotonic_on_neighbours(p, y, p) == pytest.approx(y)


def test_isotonic_clips_queries_outside_the_training_range():
    p = np.array([0.3, 0.7])
    y = np.array([0.0, 1.0])
    out = isotonic_on_neighbours(p, y, np.array([0.0, 1.0]))
    assert np.all((out >= 0.0) & (out <= 1.0))


def test_platt_recovers_an_identity_mapping_on_calibrated_input():
    rng = np.random.default_rng(1)
    p = rng.uniform(0.05, 0.95, 4000)
    y = (rng.uniform(size=p.size) < p).astype(float)
    from neighbour_calibrator import platt_on_neighbours
    a, b = platt_on_neighbours(p, y)
    # Already calibrated -> the correction should be close to a=1, b=0.
    assert a == pytest.approx(1.0, abs=0.15)
    assert b == pytest.approx(0.0, abs=0.15)
    assert apply_platt(p, a, b) == pytest.approx(p, abs=0.05)


def test_murphy_resolution_is_zero_for_any_constant_predictor():
    # This is why reliability, not ECE, is the fair calibration comparison:
    # a constant is structurally advantaged on ECE and gets no credit here.
    y = np.array([1.0] * 60 + [0.0] * 40)
    d = murphy_decomposition(np.full(100, 0.6), y)
    assert d["resolution"] == pytest.approx(0.0, abs=1e-12)
    assert d["reliability"] == pytest.approx(0.0, abs=1e-12)  # 0.6 is the base rate


def test_murphy_decomposition_reconstructs_the_brier_score():
    rng = np.random.default_rng(2)
    p = rng.uniform(size=500)
    y = (rng.uniform(size=500) < p).astype(float)
    d = murphy_decomposition(p, y, n_bins=10)
    # Brier == reliability - resolution + uncertainty, up to binning error.
    approx = d["reliability"] - d["resolution"] + d["uncertainty"]
    assert approx == pytest.approx(brier(p, y), abs=0.02)


# --- null_ece: the ECE floor a perfectly calibrated predictor still pays ------

from neighbour_calibrator import bagged_isotonic, null_ece


def test_null_ece_is_far_higher_for_a_spread_predictor_than_a_constant():
    # The whole reason raw ECE cannot be compared across these methods. Both
    # predictors below are perfectly calibrated by construction; the spread one
    # still scores several times worse purely from per-bin sampling noise.
    n = 3000
    const = np.full(n, 0.9)
    spread = np.linspace(0.05, 0.99, n)
    assert null_ece(spread, n_draws=40, seed=0) > 3 * null_ece(const, n_draws=40, seed=0)


def test_null_ece_falls_as_the_sample_grows():
    # It is a finite-sample effect, so more data must shrink it. If this ever
    # stopped holding, the floor would be measuring something else.
    spread_small = np.linspace(0.05, 0.99, 300)
    spread_big = np.linspace(0.05, 0.99, 30000)
    assert null_ece(spread_big, n_draws=20, seed=1) < null_ece(spread_small, n_draws=20, seed=1)


def test_null_ece_is_deterministic_for_a_fixed_seed():
    p = np.linspace(0.1, 0.9, 500)
    assert null_ece(p, n_draws=10, seed=7) == null_ece(p, n_draws=10, seed=7)


def test_a_genuinely_miscalibrated_predictor_sits_above_its_floor():
    # The floor test must have teeth: a predictor claiming 0.9 while being
    # right 60% of the time has to land clearly above its own floor, or the
    # test would excuse real miscalibration.
    rng = np.random.default_rng(0)
    p = np.full(3000, 0.9)
    y = (rng.uniform(size=3000) < 0.6).astype(float)
    assert ece(p, y) - null_ece(p, n_draws=40, seed=0) > 0.2


def test_bagged_isotonic_stays_monotone_and_bounded():
    rng = np.random.default_rng(3)
    p = np.sort(rng.uniform(size=250))
    y = (rng.uniform(size=250) < p).astype(float)
    out = bagged_isotonic(p, y, p, n_bags=15, seed=2)
    assert np.all(np.diff(out) >= -1e-12)
    assert np.all((out >= 0.0) & (out <= 1.0))
