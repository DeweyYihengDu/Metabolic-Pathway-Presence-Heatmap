"""Unit tests for the KEGG module-completeness evaluator."""
from mpph.modules import (
    classify_state,
    evaluate_module,
    module_completeness,
    split_steps,
    step_complete,
)


def test_split_steps_flat():
    assert split_steps("A B C") == ["A", "B", "C"]


def test_split_steps_respects_parentheses():
    assert split_steps("(A B) C") == ["(A B)", "C"]
    assert split_steps("(((K1,K2) K3),K4) K5") == ["(((K1,K2) K3),K4)", "K5"]


def test_step_or_alternatives():
    assert step_complete("K00844,K12407", {"K00844"}) is True
    assert step_complete("K00844,K12407", {"K99999"}) is False


def test_step_and_complex():
    # '+' joins essential complex subunits (AND)
    assert step_complete("K00001+K00002", {"K00001"}) is False
    assert step_complete("K00001+K00002", {"K00001", "K00002"}) is True


def test_step_optional_ignored():
    # '-' marks an optional component that must not affect completeness
    assert step_complete("K00001-K00002", {"K00001"}) is True
    assert step_complete("K00001-K00002", {"K00001", "K00002"}) is True


def test_step_nested_groups():
    defn = "(((K00134,K00150) K00927),K11389)"
    assert step_complete(defn, {"K11389"}) is True
    assert step_complete(defn, {"K00134", "K00927"}) is True
    assert step_complete(defn, {"K00134"}) is False


def test_module_completeness_fraction():
    defn = "(K00844,K12407) (K01810) K01803"
    assert module_completeness(defn, {"K00844", "K01803"}) == 2 / 3
    assert module_completeness(defn, {"K00844", "K01810", "K01803"}) == 1.0
    assert module_completeness(defn, set()) == 0.0


def test_module_completeness_empty_definition():
    assert module_completeness("", {"K00001"}) == 0.0


def test_double_dash_gap_not_counted():
    # '--' is a placeholder step (reaction with no assigned KO); it must not
    # drag completeness below 1.0 when every real step is satisfied.
    defn = "K14652 K22100 -- K01633 K13941 K22099 K00287"
    kos = {"K14652", "K22100", "K01633", "K13941", "K22099", "K00287"}
    assert module_completeness(defn, kos) == 1.0


def test_nested_module_reference_resolves():
    # A module that references sub-modules is complete iff the sub-modules are.
    module_defs = {
        "M00161": ("K00001", "cat", "Pathway"),
        "M00163": ("K00002", "cat", "Pathway"),
        "M00165": ("K00003", "cat", "Pathway"),
    }
    defn = "(M00161,M00163) M00165"
    assert module_completeness(defn, {"K00001", "K00003"}, module_defs) == 1.0
    assert module_completeness(defn, {"K00003"}, module_defs) == 0.5
    # Without module_defs the references cannot resolve -> counted absent.
    assert module_completeness(defn, {"K00001", "K00003"}) == 0.0


def test_nested_module_cycle_is_safe():
    module_defs = {"M00001": ("M00002", "c", "Pathway"),
                   "M00002": ("M00001", "c", "Pathway")}
    # Must terminate rather than recurse forever.
    assert module_completeness("M00001", {"K00001"}, module_defs) == 0.0


def test_classify_state():
    assert classify_state(1.0) == "complete"
    assert classify_state(0.5) == "partial"
    assert classify_state(0.0) == "absent"


def test_evaluate_module_evidence():
    defn = "(K00844,K12407) (K01810) K01803"
    ev = evaluate_module("M99999", defn, {"K00844", "K01803"})
    assert ev.state == "partial"
    assert ev.n_steps == 3
    assert ev.n_satisfied == 2
    assert "K00844" in ev.matched_kos
    # The unsatisfied middle step contributes its KO to missing.
    assert "K01810" in ev.missing_kos
    # Per-step detail is recorded.
    assert [s.satisfied for s in ev.steps] == [True, False, True]


def test_evaluate_module_unresolved_reference():
    ev = evaluate_module("M99999", "M00161 K00001", {"K00001"})
    assert "M00161" in ev.unresolved_references
    assert ev.parser_status == "unresolved_references"
