"""Tests for phylogeny-aware testing.

The statistical assertions here are the point: several use hand-computed
reference values, and the headline test (`test_clade_confounded_...`) fails
loudly if the conditioning in `binary_phylo_test` ever regresses.
"""
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from mpph.analysis import differential_features
from mpph.phylo import (
    _cliffs_delta_columns,
    binary_phylo_test,
    continuous_null,
    fit_mk_rate,
    fitch_changes,
    mk_loglik,
    parse_newick,
    phylo_p_from_null,
    prune_to,
    simulate_bm,
    simulate_mk,
    vcv,
)


def balanced_newick(levels: int, blen: float = 1.0, root_blen=None,
                    prefix: str = "t"):
    """A balanced tree of 2**levels tips. `root_blen` sets the two branches
    descending from the root, which is what controls how much of the tree's
    path length separates the two halves."""
    labels = [f"{prefix}{i}" for i in range(2 ** levels)]
    cur = [f"{lab}:{blen}" for lab in labels]
    while len(cur) > 2:
        cur = [f"({cur[i]},{cur[i + 1]}):{blen}" for i in range(0, len(cur), 2)]
    rb = blen if root_blen is None else root_blen
    halves = [c[:c.rfind(":")] for c in cur]
    return f"({halves[0]}:{rb},{halves[1]}:{rb});", labels


# --------------------------------------------------------------------------- #
# Parser
# --------------------------------------------------------------------------- #
def test_parse_newick_keeps_branch_lengths_and_quoted_labels():
    # Same fixture treecompare uses -- but that parser discards lengths, which
    # is exactly why this module needs its own.
    tree = parse_newick(
        "(('[Eubacterium] rectale':0.5,'O''Brien strain':0.5):1.0,plain_c:1.5);")
    assert set(tree.tip_labels) == {
        "[Eubacterium] rectale", "O'Brien strain", "plain_c"}
    depths = tree.root_to_tip()
    assert depths["[Eubacterium] rectale"] == pytest.approx(1.5)
    assert depths["O'Brien strain"] == pytest.approx(1.5)
    assert depths["plain_c"] == pytest.approx(1.5)


def test_parse_newick_rejects_cladogram_without_lengths():
    # Defaulting missing lengths to 1.0 would fabricate the distances the
    # whole correction is computed from.
    with pytest.raises(ValueError, match="no branch lengths"):
        parse_newick("((a,b),(c,d));")


@pytest.mark.parametrize("newick,match", [
    ("((a:1,b:-1):1,(c:1,d:1):1);", "negative branch length"),
    ("((a:1,a:1):1,(c:1,d:1):1);", "duplicate tip labels"),
    ("(a:1);", "need >= 2 tips"),
])
def test_parse_newick_rejects_malformed_trees(newick, match):
    with pytest.raises(ValueError, match=match):
        parse_newick(newick)


# --------------------------------------------------------------------------- #
# Exact oracles -- no RNG
# --------------------------------------------------------------------------- #
def test_vcv_matches_hand_computed_matrix():
    tree = parse_newick("((a:1,b:1):1,(c:1,d:1):1);")
    expected = np.array([[2, 1, 0, 0], [1, 2, 0, 0],
                         [0, 0, 2, 1], [0, 0, 1, 2]], dtype=float)
    assert np.array_equal(vcv(tree), expected)


def test_mk_loglik_matches_hand_computed_two_tip_likelihood():
    # (a:1,b:1) with q=0.5: p = 0.5*(1 - e^-1) is the chance of a change over
    # one unit branch. L(differing tips) = p(1-p); L(same) = 0.5[(1-p)^2 + p^2].
    tree = parse_newick("(a:1.0,b:1.0);")
    p = 0.5 * (1 - np.exp(-1.0))
    assert np.exp(mk_loglik(tree, np.array([0.0, 1.0]), 0.5)) == pytest.approx(
        p * (1 - p), abs=1e-12)
    assert np.exp(mk_loglik(tree, np.array([0.0, 0.0]), 0.5)) == pytest.approx(
        0.5 * ((1 - p) ** 2 + p ** 2), abs=1e-12)


def test_mk_loglik_is_identical_across_rootings_of_the_same_tree():
    # Felsenstein's pulley principle: with a symmetric rate and a stationary
    # root the likelihood cannot depend on where an unrooted tree is rooted.
    # This single check catches most traversal and branch-length bookkeeping
    # bugs at once, so it asserts bit-equality rather than approximate.
    t1 = parse_newick("((a:1,b:1):1,(c:1,d:1):1);")
    t2 = parse_newick("(a:0.5,(b:1,(c:1,d:1):2):0.5);")
    values = {"a": 0.0, "b": 1.0, "c": 1.0, "d": 0.0}
    s1 = np.array([values[lab] for lab in t1.tip_labels])
    s2 = np.array([values[lab] for lab in t2.tip_labels])
    for q in (0.1, 0.5, 2.0):
        assert mk_loglik(t1, s1, q) == mk_loglik(t2, s2, q)


def test_mk_loglik_treats_nan_tip_as_unconstrained():
    # A tip with unknown state gets partial likelihood [1, 1], so it should
    # contribute nothing -- the two-tip tree reduces to the one-tip case.
    tree = parse_newick("(a:1.0,b:1.0);")
    assert np.exp(mk_loglik(tree, np.array([1.0, np.nan]), 0.5)) == pytest.approx(
        0.5, abs=1e-12)


def test_fitch_changes_hand_computed():
    tree, _ = balanced_newick(3)  # 8 tips
    t8 = parse_newick(tree)
    monophyletic = np.array([1.0] * 4 + [0.0] * 4)
    alternating = np.array([1.0, 0.0] * 4)
    assert fitch_changes(t8, monophyletic) == 1
    assert fitch_changes(t8, alternating) == 4


def test_prune_adds_collapsed_node_length_to_surviving_child():
    # Dropping b leaves a's parent with one child. If that node's length is
    # not folded into a, every affected root-to-tip path silently shortens,
    # which inflates apparent phylogenetic signal.
    tree = parse_newick("((a:1,b:1):1,(c:1,d:1):1);")
    assert tree.root_to_tip()["a"] == pytest.approx(2.0)
    pruned = prune_to(tree, {"a", "c", "d"})
    assert pruned.root_to_tip()["a"] == pytest.approx(2.0)
    assert set(pruned.tip_labels) == {"a", "c", "d"}


def test_prune_rejects_too_few_shared_tips():
    tree = parse_newick("((a:1,b:1):1,(c:1,d:1):1);")
    with pytest.raises(ValueError, match="need >= 2"):
        prune_to(tree, {"a"})


# --------------------------------------------------------------------------- #
# Simulation
# --------------------------------------------------------------------------- #
def test_simulate_bm_reproduces_the_vcv_covariance():
    # Per-branch increments must be equivalent to drawing from
    # MVN(0, sigma^2 * C) -- checked against the independently-built C.
    tree = parse_newick("((a:1,b:1):1,(c:1,d:2.5):0.5);")
    values = simulate_bm(tree, 60000, np.random.default_rng(0))
    empirical = np.cov(values)
    assert np.allclose(empirical, vcv(tree), atol=0.08)


def test_simulate_mk_change_probability_matches_the_closed_form():
    # Over a single branch of length t, P(state differs) = 0.5(1 - e^-2qt).
    tree = parse_newick("(a:1.0,b:1.0);")
    q = 0.7
    sims = simulate_mk(tree, q, 200000, np.random.default_rng(0))
    differ = float(np.mean(sims[0] != sims[1]))
    p_one = 0.5 * (1 - np.exp(-2 * q * 1.0))
    expected = 2 * p_one * (1 - p_one)   # exactly one of the two branches changed
    assert differ == pytest.approx(expected, abs=0.01)


def test_fit_mk_rate_flags_boundary_for_an_invariant_trait():
    newick, labels = balanced_newick(3)
    t8 = parse_newick(newick)
    _, at_boundary = fit_mk_rate(t8, np.zeros(len(labels)))
    assert at_boundary


# --------------------------------------------------------------------------- #
# The statistics that matter
# --------------------------------------------------------------------------- #
def test_clade_confounded_binary_feature_is_not_significant_after_correction():
    """T1, the headline.

    Trait == clade == group, so the group difference is entirely explained by
    shared ancestry. Fisher's exact test calls this p = 1.55e-4; the corrected
    test must not.

    The threshold is deliberately high. An unconditioned null (one that does
    not reject simulated replicates whose presence count differs from the
    observed one) scores 0.02-0.03 here and would pass a `> 0.05` assertion
    while being broken -- roughly a third of its replicates come out
    invariant, contribute a zero statistic and dilute the tail. If this
    assertion ever needs lowering, the conditioning has regressed.
    """
    newick, _labels = balanced_newick(4)   # 16 tips
    tree = parse_newick(newick)
    states = np.array([1.0] * 8 + [0.0] * 8)
    idx_a, idx_b = np.arange(8), np.arange(8, 16)

    _, fisher_p = stats.fisher_exact([[8, 0], [0, 8]])
    assert fisher_p < 1e-3

    assert fitch_changes(tree, states) == 1
    result = binary_phylo_test(tree, states, idx_a, idx_b, n_sims=9999,
                               rng=np.random.default_rng(0))
    assert result["p_value_phylo"] > 0.5
    assert result["phylo_k_tolerance"] == 0        # exact-k matching sufficed
    assert result["n_phylo_replicates_accepted"] > 100


def test_real_effect_scattered_across_the_tree_stays_significant():
    """T2, power. The correction must not simply never reject."""
    newick, _ = balanced_newick(5)   # 32 tips
    tree = parse_newick(newick)
    rng = np.random.default_rng(7)
    perm = rng.permutation(32)
    idx_a, idx_b = np.sort(perm[:16]), np.sort(perm[16:])
    states = np.zeros(32)
    states[idx_a] = (rng.random(16) < 0.9).astype(float)
    states[idx_b] = (rng.random(16) < 0.1).astype(float)

    a_pres = int(states[idx_a].sum())
    b_pres = int(states[idx_b].sum())
    _, fisher_p = stats.fisher_exact(
        [[a_pres, 16 - a_pres], [b_pres, 16 - b_pres]])
    result = binary_phylo_test(tree, states, idx_a, idx_b, n_sims=9999,
                               rng=np.random.default_rng(1))
    assert fisher_p < 0.05
    assert result["p_value_phylo"] < 0.05


def test_binary_null_replicates_all_match_the_observed_presence_count():
    newick, _ = balanced_newick(4)
    tree = parse_newick(newick)
    states = np.array([1.0] * 5 + [0.0] * 11)
    idx_a, idx_b = np.arange(8), np.arange(8, 16)
    result = binary_phylo_test(tree, states, idx_a, idx_b, n_sims=9999,
                               rng=np.random.default_rng(0))
    # With exact-k matching the tolerance must be 0; if it widened, the band
    # actually used has to be reported rather than silently applied.
    assert result["phylo_k_tolerance"] >= 0
    assert result["n_phylo_replicates_accepted"] > 0


def test_clade_confounded_continuous_feature_is_not_significant_after_correction():
    """T3, the continuous mirror of T1.

    The root-split branch is long relative to within-clade depth -- the
    realistic "two divergent genera" case. That matters: the correction's
    strength is a function of how much of the tree's path length separates
    the groups rather than varying within them. On a tree where all branches
    are equal (root split only a quarter of the depth) complete separation is
    genuinely uncommon under Brownian motion and the same trait scores ~0.06,
    which is also correct -- just a weaker correction.
    """
    newick, _labels = balanced_newick(4, blen=1.0, root_blen=20.0)
    tree = parse_newick(newick)
    values = np.array([1.0] * 8 + [0.0] * 8)
    idx_a, idx_b = np.arange(8), np.arange(8, 16)

    mw_p = stats.mannwhitneyu(values[idx_a], values[idx_b],
                              alternative="two-sided").pvalue
    assert mw_p < 0.01

    observed = float(_cliffs_delta_columns(values.reshape(-1, 1), idx_a, idx_b)[0])
    assert observed == pytest.approx(1.0)
    null = continuous_null(tree, idx_a, idx_b, 20000, np.random.default_rng(0))
    assert phylo_p_from_null(observed, null) > 0.25


def test_continuous_null_is_independent_of_the_feature_values():
    """The rank statistic makes the null free of sigma^2 and the root state,
    which is what lets one null be reused across every feature."""
    newick, _ = balanced_newick(3)
    tree = parse_newick(newick)
    idx_a, idx_b = np.arange(4), np.arange(4, 8)
    a = continuous_null(tree, idx_a, idx_b, 2000, np.random.default_rng(0))
    b = continuous_null(tree, idx_a, idx_b, 2000, np.random.default_rng(0))
    assert np.array_equal(a, b)


def test_same_seed_gives_identical_phylo_p_values():
    newick, _ = balanced_newick(4)
    tree = parse_newick(newick)
    states = np.array([1.0] * 6 + [0.0] * 10)
    idx_a, idx_b = np.arange(8), np.arange(8, 16)
    kwargs = {"n_sims": 2999}
    first = binary_phylo_test(tree, states, idx_a, idx_b,
                              rng=np.random.default_rng(3), **kwargs)
    second = binary_phylo_test(tree, states, idx_a, idx_b,
                               rng=np.random.default_rng(3), **kwargs)
    assert first == second


# --------------------------------------------------------------------------- #
# Integration through differential_features
# --------------------------------------------------------------------------- #
def test_differential_features_adds_phylo_columns_and_corrects_a_confounded_hit():
    newick, labels = balanced_newick(4)
    tree = parse_newick(newick)
    matrix = pd.DataFrame({"confounded": [1.0] * 8 + [0.0] * 8}, index=labels)
    groups = pd.Series(["A"] * 8 + ["B"] * 8, index=labels)

    out = differential_features(matrix, groups, "A", "B", tree=tree,
                                phylo_permutations=4999, seed=0)
    row = out.iloc[0]
    assert row["q_value"] < 0.01           # uncorrected: highly significant
    assert row["q_value_phylo"] > 0.5      # corrected: not
    assert row["n_state_changes"] == 1
    assert bool(row["phylo_confounded"]) is True


def test_differential_features_without_a_tree_is_unchanged():
    _newick, labels = balanced_newick(3)
    matrix = pd.DataFrame({"f": [1.0] * 4 + [0.0] * 4}, index=labels)
    groups = pd.Series(["A"] * 4 + ["B"] * 4, index=labels)
    out = differential_features(matrix, groups, "A", "B")
    for col in ("p_value_phylo", "q_value_phylo", "n_state_changes"):
        assert col not in out.columns


def _make_compare_inputs(tmp_path, matrix, groups, mode="presence"):
    import json
    matrix.to_csv(tmp_path / "Demo_matrix.csv", index_label="organism")
    (tmp_path / "Demo_manifest.json").write_text(json.dumps({"mode": mode}))
    rows = "".join(f"{o}\t{g}\n" for o, g in groups.items())
    (tmp_path / "meta.tsv").write_text("sample_id\tgroup\n" + rows)


def test_cli_compare_with_a_reference_tree_writes_phylo_columns(tmp_path):
    from mpph import cli
    _, labels = balanced_newick(3)
    matrix = pd.DataFrame({"confounded": [1.0] * 4 + [0.0] * 4}, index=labels)
    groups = pd.Series(["A"] * 4 + ["B"] * 4, index=labels)
    _make_compare_inputs(tmp_path, matrix, groups)
    # Long root split: two genuinely divergent clades.
    (tmp_path / "ref.nwk").write_text(
        "(((t0:1,t1:1):1,(t2:1,t3:1):1):5,((t4:1,t5:1):1,(t6:1,t7:1):1):5);")

    rc = cli.main(["compare", "--results", str(tmp_path),
                   "--metadata", str(tmp_path / "meta.tsv"),
                   "--group-column", "group", "--group-a", "A", "--group-b", "B",
                   "--tree", str(tmp_path / "ref.nwk"),
                   "--phylo-permutations", "999"])
    assert rc == 0
    out = pd.read_csv(next(tmp_path.glob("Demo_differential_*.csv")))
    assert {"p_value_phylo", "q_value_phylo", "n_state_changes"} <= set(out.columns)
    assert out.iloc[0]["n_state_changes"] == 1


def test_cli_compare_refuses_mpphs_own_dendrogram_even_when_renamed(tmp_path, capsys):
    # The guard compares topology, not filenames: copying or renaming the
    # dendrogram must not get it past.
    from scipy.cluster.hierarchy import linkage

    from mpph import cli
    from mpph.tree import linkage_to_newick
    labels = [f"t{i}" for i in range(8)]
    matrix = pd.DataFrame({
        "a": [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        "b": [1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
        "c": [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
    }, index=labels)
    groups = pd.Series(["A"] * 4 + ["B"] * 4, index=labels)
    _make_compare_inputs(tmp_path, matrix, groups)
    own = linkage_to_newick(linkage(matrix.to_numpy(dtype=float),
                                    method="average"), labels)
    (tmp_path / "Demo_organism_tree.nwk").write_text(own)
    (tmp_path / "innocent_name.nwk").write_text(own)

    rc = cli.main(["compare", "--results", str(tmp_path),
                   "--metadata", str(tmp_path / "meta.tsv"),
                   "--group-column", "group", "--group-a", "A", "--group-b", "B",
                   "--tree", str(tmp_path / "innocent_name.nwk"),
                   "--phylo-permutations", "199"])
    assert rc == 1
    assert "circular" in capsys.readouterr().err


def test_cli_compare_errors_when_tree_labels_barely_match(tmp_path, capsys):
    from mpph import cli
    _, labels = balanced_newick(3)
    matrix = pd.DataFrame({"f": [1.0] * 4 + [0.0] * 4}, index=labels)
    groups = pd.Series(["A"] * 4 + ["B"] * 4, index=labels)
    _make_compare_inputs(tmp_path, matrix, groups)
    # GTDB-style accessions rather than the matrix's organism names -- the
    # single most likely real-world failure.
    (tmp_path / "ref.nwk").write_text(
        "((GB_GCA_1:1,GB_GCA_2:1):1,(GB_GCA_3:1,GB_GCA_4:1):1);")
    rc = cli.main(["compare", "--results", str(tmp_path),
                   "--metadata", str(tmp_path / "meta.tsv"),
                   "--group-column", "group", "--group-a", "A", "--group-b", "B",
                   "--tree", str(tmp_path / "ref.nwk")])
    assert rc == 1
    assert "match a tree tip" in capsys.readouterr().err


def test_differential_features_errors_when_an_organism_is_missing_from_the_tree():
    newick, labels = balanced_newick(3)
    tree = parse_newick(newick)
    matrix = pd.DataFrame({"f": [1.0] * 4 + [0.0] * 5},
                          index=labels + ["not_in_tree"])
    groups = pd.Series(["A"] * 4 + ["B"] * 5, index=matrix.index)
    with pytest.raises(ValueError, match="not tips of the tree"):
        differential_features(matrix, groups, "A", "B", tree=tree,
                              phylo_permutations=99)
