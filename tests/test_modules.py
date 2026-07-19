"""Unit tests for the KEGG module-completeness evaluator."""
import math

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
    # Without module_defs the references can't be resolved either way -> every
    # step is undetermined, not "confirmed absent" (NaN, never a fabricated 0.0).
    assert math.isnan(module_completeness(defn, {"K00001", "K00003"}))


def test_nested_module_cycle_is_unknown_not_absent():
    module_defs = {"M00001": ("M00002", "c", "Pathway"),
                   "M00002": ("M00001", "c", "Pathway")}
    # Must terminate rather than recurse forever, and a pure cycle carries no
    # real evidence either way -> NaN ("unknown"), not a fabricated 0.0.
    assert math.isnan(module_completeness("M00001", {"K00001"}, module_defs))


def test_step_or_with_known_hit_ignores_unresolved_reference():
    # The OR is already satisfied by a present KO -- the unresolved module
    # reference can't change that outcome, so this must be determined True,
    # not "unknown", even though module_defs doesn't cover M00161.
    assert step_complete("K00001,M00161", {"K00001"}) is True


def test_step_and_with_known_miss_ignores_unresolved_reference():
    # The AND is already false (the KO is missing) regardless of the
    # unresolved reference -- false dominates AND, so this is determined
    # False, not "unknown".
    assert step_complete("K00001+M00161", {"K99999"}) is False


def test_step_genuinely_depends_on_unresolved_reference_is_unknown():
    # Whether this step is satisfied genuinely hinges on M00161, which we
    # have no definition for -- must be unknown (None), never coerced False.
    assert step_complete("K00001+M00161", {"K00001"}) is None
    assert step_complete("K00001,M00161", {"K99999"}) is None


def test_module_completeness_excludes_unknown_steps_from_ratio():
    # 2 determinable steps (1 satisfied), 1 genuinely-unknown step (an
    # unresolved reference this OR can't route around) -- the ratio must be
    # computed over the 2 determinable steps only (1/2), not 1/3.
    defn = "K00001 K99999 K00002,M00161"
    assert module_completeness(defn, {"K00001"}) == 0.5


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
    # The AND step genuinely depends on the unresolved M00161 -> unknown.
    assert ev.steps[0].satisfied is None


def test_evaluate_module_state_unknown_when_wholly_undetermined():
    # Every step is a bare, unresolved module reference -- nothing is
    # determinable, so the module's state must be "unknown", never a
    # fabricated "absent" derived from treating unresolved as failed.
    ev = evaluate_module("M99999", "M00161 M00162", {"K00001"})
    assert math.isnan(ev.score)
    assert ev.state == "unknown"
    assert ev.n_steps == 0
    assert ev.n_satisfied == 0


def test_evaluate_module_unresolved_reference_not_reported_when_outcome_determined():
    # The OR is already satisfied by a present KO, so M00161 never actually
    # influenced the outcome -- it should not be reported as "unresolved".
    ev = evaluate_module("M99999", "K00001,M00161", {"K00001"})
    assert ev.steps[0].satisfied is True
    assert ev.unresolved_references == set()
    assert ev.parser_status == "valid"
