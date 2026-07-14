"""Unit tests for the KEGG module-completeness evaluator."""
from mpph.modules import module_completeness, split_steps, step_complete


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
