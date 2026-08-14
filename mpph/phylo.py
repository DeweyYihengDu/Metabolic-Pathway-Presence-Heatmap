"""Phylogeny-aware machinery for `mpph compare`: a Newick parser that keeps
branch lengths, trait simulation under Brownian motion and a symmetric Mk
model, and a parsimony diagnostic.

Why this exists: `differential_features` treats every genome as an
independent observation, but close relatives share features by descent, so
its p-values are anti-conservative. Measured on a 32-tip balanced tree with
the two groups being the two root-split clades and traits simulated with *no*
group effect, Fisher's exact test rejects at 27.7% for a nominal 5% test.
The fix is a null distribution generated on the tree instead of one that
assumes exchangeable samples.

Two properties of this design are worth stating up front because they are not
obvious and they are what keep it honest:

* **Root-invariance.** With a symmetric rate and a stationary root the Mk
  likelihood is identical whichever way an unrooted tree is rooted
  (Felsenstein's pulley principle), and every Brownian contrast variance
  ``C_ii + C_jj - 2*C_ij`` is likewise invariant. Rank and difference
  statistics see the data only through contrasts, so no rerooting step is
  needed for either branch.
* **The continuous null has no fitted parameters.** Under Brownian motion a
  trait is ``sigma * Z + root_state``; a *rank* statistic is invariant to
  positive scaling and to translation, so its null distribution depends on
  neither. That removes an optimizer, a parameter-estimation bias argument,
  and a whole class of misspecification failure from that branch.

``mpph``'s own dendrogram must never be used as the tree here -- it is built
from the very features being tested, so correcting those tests with it is
circular. The caller is responsible for enforcing that; see
:func:`mpph.analysis.differential_features`.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

# Newick quoting: a quoted label ('...' with '' as an escaped internal quote),
# optionally followed by a :branch_length, is one token. mpph.tree's writer
# quotes labels rather than stripping characters from them, so any parser here
# has to round-trip that. Defined here and imported by treecompare so the two
# parsers cannot drift apart on what counts as a label.
_QUOTED_TOKEN = re.compile(r"'(?:[^']|'')*'(?::[^(),;]*)?")
_QUOTED_NAME = re.compile(r"^'((?:[^']|'')*)'")
_TOKEN = re.compile(_QUOTED_TOKEN.pattern + r"|[(),;]|[^(),;]+")


@dataclass(frozen=True)
class PhyloTree:
    """A rooted tree with branch lengths, stored as parallel arrays.

    Nodes are ordered so that every parent has a lower index than all of its
    children, which lets the root-to-tip simulation passes run as a single
    forward loop over ``range(1, n_nodes)`` -- no recursion, and no risk of
    hitting Python's recursion limit on a deep tree.

    ``parent[0] == -1`` marks the root; ``blen[0]`` is 0 by convention (a root
    branch has no meaning here and is ignored by every method).
    """
    parent: np.ndarray          # int, parent index per node (-1 at the root)
    blen: np.ndarray            # float, branch length above each node
    is_tip: np.ndarray          # bool
    labels: list                # str per node ("" for unlabelled internals)
    children: list              # list[list[int]]

    @property
    def n_nodes(self) -> int:
        return len(self.parent)

    @property
    def tip_indices(self) -> np.ndarray:
        return np.flatnonzero(self.is_tip)

    @property
    def tip_labels(self) -> list:
        return [self.labels[i] for i in self.tip_indices]

    def tip_index(self) -> dict:
        """``{tip label: node index}``."""
        return {self.labels[i]: int(i) for i in self.tip_indices}

    def root_to_tip(self) -> dict:
        """``{tip label: summed branch length from the root}`` -- the check
        that catches a botched degree-2 collapse in :func:`prune_to`."""
        depth = np.zeros(self.n_nodes)
        for i in range(1, self.n_nodes):
            depth[i] = depth[self.parent[i]] + self.blen[i]
        return {self.labels[i]: float(depth[i]) for i in self.tip_indices}


def parse_newick(text: str) -> PhyloTree:
    """Parse a Newick string, **keeping branch lengths**.

    ``mpph.treecompare`` also parses Newick but deliberately keeps only
    topology (it computes Robinson-Foulds on bipartitions), so it cannot be
    reused here: every method in this module needs branch lengths.

    A tree without branch lengths (a cladogram) is rejected rather than
    defaulted to 1.0. ``mpph``'s own writer always emits lengths, so a
    length-free tree is necessarily a user's file, and silently inventing
    lengths would produce confident nonsense from it.
    """
    tokens = [t for t in _TOKEN.findall(text.strip()) if t.strip()]
    if not tokens:
        raise ValueError("empty Newick string")

    parent: list = []
    blen: list = []
    labels: list = []
    children: list = []
    stack: list = []
    last_closed = None
    seen_length = False

    def _add_node(par: int) -> int:
        idx = len(parent)
        parent.append(par)
        blen.append(0.0)
        labels.append("")
        children.append([])
        if par >= 0:
            children[par].append(idx)
        return idx

    def _split_label(tok: str):
        """``name`` / ``name:len`` / ``'quoted name':len`` -> (name, length)."""
        quoted = _QUOTED_NAME.match(tok)
        if quoted:
            name = quoted.group(1).replace("''", "'")
            rest = tok[quoted.end():]
            length = rest[1:].strip() if rest.startswith(":") else None
        elif ":" in tok:
            name, _, raw = tok.partition(":")
            name, length = name.strip(), raw.strip()
        else:
            name, length = tok.strip(), None
        return name, length

    for tok in tokens:
        if tok == "(":
            # The outermost '(' creates the root (parent -1); nested ones
            # create internal children.
            node = _add_node(stack[-1] if stack else -1)
            stack.append(node)
            last_closed = None
        elif tok == ",":
            last_closed = None
        elif tok == ")":
            if not stack:
                raise ValueError("malformed Newick: unbalanced ')'")
            last_closed = stack.pop()
        elif tok == ";":
            break
        else:
            name, length = _split_label(tok)
            if last_closed is not None:
                # A label and/or :length attached to the node just closed.
                node = last_closed
                last_closed = None
            elif stack:
                node = _add_node(stack[-1])   # a leaf
            else:
                raise ValueError(f"malformed Newick: stray token {tok!r}")
            if name:
                labels[node] = name
            if length is not None:
                try:
                    blen[node] = float(length)
                except ValueError as exc:
                    raise ValueError(
                        f"malformed branch length {length!r}") from exc
                seen_length = True

    if stack:
        raise ValueError("malformed Newick: unbalanced parentheses")
    if not parent:
        raise ValueError("malformed Newick: no tree found")

    parent_arr = np.asarray(parent, dtype=int)
    blen_arr = np.asarray(blen, dtype=float)
    child_counts = np.asarray([len(c) for c in children], dtype=int)
    is_tip = child_counts == 0

    if not seen_length:
        raise ValueError(
            "this Newick tree has no branch lengths (a cladogram). Branch "
            "lengths are required: defaulting them to 1.0 would silently "
            "fabricate the evolutionary distances the correction depends on.")
    if np.any(blen_arr[1:] < 0):
        raise ValueError("negative branch length in Newick tree")
    tip_labels = [labels[i] for i in np.flatnonzero(is_tip)]
    if any(not lab for lab in tip_labels):
        raise ValueError("unlabelled tip in Newick tree")
    if len(set(tip_labels)) != len(tip_labels):
        dupes = sorted({lab for lab in tip_labels if tip_labels.count(lab) > 1})
        raise ValueError(f"duplicate tip labels in Newick tree: {dupes[:5]}")
    if len(tip_labels) < 2:
        raise ValueError("need >= 2 tips")

    tree = PhyloTree(parent=parent_arr, blen=blen_arr, is_tip=is_tip,
                     labels=labels, children=children)
    return _reorder_parents_first(tree)


def _reorder_parents_first(tree: PhyloTree) -> PhyloTree:
    """Renumber nodes so every parent precedes its children."""
    order: list = []
    stack = [0]
    while stack:
        node = stack.pop()
        order.append(node)
        stack.extend(reversed(tree.children[node]))
    remap = {old: new for new, old in enumerate(order)}
    parent = np.asarray(
        [-1 if tree.parent[o] < 0 else remap[int(tree.parent[o])] for o in order],
        dtype=int)
    blen = np.asarray([tree.blen[o] for o in order], dtype=float)
    is_tip = np.asarray([tree.is_tip[o] for o in order], dtype=bool)
    labels = [tree.labels[o] for o in order]
    children = [[remap[c] for c in tree.children[o]] for o in order]
    return PhyloTree(parent=parent, blen=blen, is_tip=is_tip, labels=labels,
                     children=children)


def prune_to(tree: PhyloTree, keep) -> PhyloTree:
    """Restrict ``tree`` to the tips in ``keep``.

    When dropping tips leaves an internal node with a single child, that node
    is collapsed and **its branch length is added to the surviving child**.
    Skipping that addition is the classic silent-wrong-answer bug here: it
    shortens every affected root-to-tip path, which inflates apparent
    phylogenetic signal and makes the correction quietly too conservative.
    :meth:`PhyloTree.root_to_tip` before/after is the assertion that catches it.
    """
    keep = {str(k) for k in keep}
    present = [lab for lab in tree.tip_labels if lab in keep]
    if len(present) < 2:
        raise ValueError(
            f"pruning left {len(present)} tip(s); need >= 2 shared with the tree")

    # Bottom-up: which nodes still subtend a kept tip?
    alive = np.zeros(tree.n_nodes, dtype=bool)
    for i in range(tree.n_nodes - 1, -1, -1):
        if tree.is_tip[i]:
            alive[i] = tree.labels[i] in keep
        else:
            alive[i] = any(alive[c] for c in tree.children[i])

    parent: list = []
    blen: list = []
    labels: list = []
    children: list = []
    index_of: dict = {}

    def _emit(node: int, accumulated: float, new_parent: int) -> None:
        """Copy `node` into the new tree, folding in `accumulated` length from
        any degree-2 ancestors that were collapsed on the way down."""
        live_children = [c for c in tree.children[node] if alive[c]]
        length = accumulated + (tree.blen[node] if new_parent >= 0 else 0.0)
        if len(live_children) == 1 and new_parent >= 0:
            # Degree-2: don't emit this node, carry its length to the child.
            _emit(live_children[0], length, new_parent)
            return
        idx = len(parent)
        index_of[node] = idx
        parent.append(new_parent)
        blen.append(length)
        labels.append(tree.labels[node])
        children.append([])
        if new_parent >= 0:
            children[new_parent].append(idx)
        for child in live_children:
            _emit(child, 0.0, idx)

    root = 0
    # A root with one live child is itself degree-2: descend past it.
    while not tree.is_tip[root]:
        live = [c for c in tree.children[root] if alive[c]]
        if len(live) != 1:
            break
        root = live[0]
    _emit(root, 0.0, -1)

    parent_arr = np.asarray(parent, dtype=int)
    blen_arr = np.asarray(blen, dtype=float)
    blen_arr[0] = 0.0
    is_tip = np.asarray([len(c) == 0 for c in children], dtype=bool)
    return _reorder_parents_first(PhyloTree(
        parent=parent_arr, blen=blen_arr, is_tip=is_tip, labels=labels,
        children=children))


def vcv(tree: PhyloTree) -> np.ndarray:
    """Brownian covariance between tips: shared root-to-MRCA path length.

    Only used as a test oracle and for documentation -- the simulation paths
    below never form this matrix (they add per-branch increments instead,
    which is O(n) rather than O(n^2) and needs no decomposition).
    """
    tips = tree.tip_indices
    depth = np.zeros(tree.n_nodes)
    for i in range(1, tree.n_nodes):
        depth[i] = depth[tree.parent[i]] + tree.blen[i]

    ancestors = []
    for t in tips:
        chain = set()
        node = int(t)
        while node >= 0:
            chain.add(node)
            node = int(tree.parent[node])
        ancestors.append(chain)

    n = len(tips)
    out = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            shared = ancestors[i] & ancestors[j]
            mrca = max(shared, key=lambda node: depth[node])
            out[i, j] = out[j, i] = depth[mrca]
    return out


def simulate_bm(tree: PhyloTree, n_sims: int, rng) -> np.ndarray:
    """Simulate Brownian motion down the tree; returns ``(n_tips, n_sims)``.

    The rate is fixed at 1 and the root at 0 **on purpose**: the caller scores
    these with a rank statistic, which is invariant to positive scaling and to
    translation, so the null distribution does not depend on either. That is
    what lets the continuous branch avoid fitting any parameter at all.

    Increments are added per branch rather than drawing from
    ``MVN(0, sigma^2 * vcv(tree))`` -- the two are equivalent (verified against
    :func:`vcv` in the tests) but this is O(nodes) with no decomposition.
    """
    values = np.zeros((tree.n_nodes, n_sims))
    sd = np.sqrt(tree.blen)
    for i in range(1, tree.n_nodes):
        values[i] = values[tree.parent[i]]
        if sd[i] > 0:
            values[i] += rng.normal(0.0, sd[i], n_sims)
    return values[tree.tip_indices]


def simulate_mk(tree: PhyloTree, q: float, n_sims: int, rng) -> np.ndarray:
    """Simulate a symmetric 2-state Mk trait; returns ``(n_tips, n_sims)`` uint8.

    Over a branch of length ``t`` the chance of ending in the other state is
    ``0.5 * (1 - exp(-2*q*t))`` (the eigenvalues of ``q*[[-1,1],[1,-1]]`` are
    0 and -2q). The root is drawn from the stationary distribution [0.5, 0.5].
    """
    p_change = 0.5 * (1.0 - np.exp(-2.0 * q * tree.blen))
    states = np.empty((tree.n_nodes, n_sims), dtype=np.uint8)
    states[0] = (rng.random(n_sims) < 0.5).astype(np.uint8)
    for i in range(1, tree.n_nodes):
        flip = (rng.random(n_sims) < p_change[i]).astype(np.uint8)
        states[i] = states[tree.parent[i]] ^ flip
    return states[tree.tip_indices]


def mk_loglik(tree: PhyloTree, states: np.ndarray, q: float) -> float:
    """Log-likelihood of ``states`` under a symmetric 2-state Mk model.

    Felsenstein's pruning algorithm. Each node's partial-likelihood vector is
    normalised and the log of the scale factor accumulated -- a plain product
    underflows float64 well before a few hundred tips.

    A tip whose state is NaN gets partial likelihood ``[1, 1]``, which is the
    mathematically correct treatment of "unknown" (it contributes no
    constraint) and means missing data needs no per-feature tree pruning.
    """
    if q <= 0:
        raise ValueError("Mk rate q must be positive")
    partials = np.ones((tree.n_nodes, 2))
    for pos, node in enumerate(tree.tip_indices):
        value = states[pos]
        if not np.isnan(value):
            state = int(value > 0.5)
            partials[node] = [0.0, 0.0]
            partials[node][state] = 1.0

    exp_term = np.exp(-2.0 * q * tree.blen)
    p_same = 0.5 + 0.5 * exp_term
    p_diff = 0.5 - 0.5 * exp_term

    log_scale = 0.0
    for i in range(tree.n_nodes - 1, -1, -1):
        if tree.is_tip[i]:
            continue
        acc = np.ones(2)
        for child in tree.children[i]:
            cp = partials[child]
            acc *= [p_same[child] * cp[0] + p_diff[child] * cp[1],
                    p_diff[child] * cp[0] + p_same[child] * cp[1]]
        total = acc.sum()
        if total <= 0:
            return -np.inf
        partials[i] = acc / total
        log_scale += np.log(total)
    return float(np.log(0.5 * partials[0].sum()) + log_scale)


def fit_mk_rate(tree: PhyloTree, states: np.ndarray,
                bounds=(1e-4, 50.0)) -> tuple:
    """ML estimate of the symmetric Mk rate; returns ``(q_hat, at_boundary)``.

    Optimised over ``log q`` within ``bounds``. The rate pins to the lower
    bound for a near-invariant trait and can pin to the upper bound for a
    maximally scattered one; ``at_boundary`` says so rather than letting a
    null distribution be built silently from a bound.

    Fitting a single rate with a stationary root and no group term *is*
    fitting the null model, which makes the simulation that uses it a
    parametric bootstrap with a plug-in nuisance parameter. Both directions
    of bias are safe here: a real clade-confounded effect drives q down and a
    real scattered effect drives it up, and both widen the null.
    """
    from scipy.optimize import minimize_scalar

    lo, hi = np.log(bounds[0]), np.log(bounds[1])
    result = minimize_scalar(
        lambda log_q: -mk_loglik(tree, states, float(np.exp(log_q))),
        bounds=(lo, hi), method="bounded")
    q_hat = float(np.exp(result.x))
    # Relative, not absolute: `bounded` stops within its own xatol of the
    # bound rather than exactly on it, so an absolute tolerance tighter than
    # that convergence tolerance never fires and the flag silently never
    # trips (it lands at e.g. 1.00000061e-4 against a 1e-4 bound).
    at_boundary = bool(q_hat <= bounds[0] * 1.01 or q_hat >= bounds[1] * 0.99)
    return q_hat, at_boundary


def _cliffs_delta_columns(values: np.ndarray, idx_a: np.ndarray,
                          idx_b: np.ndarray) -> np.ndarray:
    """Cliff's delta per column of a ``(n_tips, n_sims)`` array.

    Computed from ranks (``delta = 2U/(n*m) - 1``) so it vectorises across
    simulations; simulated Brownian values are continuous, so exact ties have
    probability zero and mid-ranks are unnecessary.
    """
    a, b = values[idx_a], values[idx_b]
    n_a, n_b = len(idx_a), len(idx_b)
    combined = np.concatenate([a, b], axis=0)
    order = np.argsort(combined, axis=0, kind="stable")
    ranks = np.empty_like(order)
    rows = np.arange(combined.shape[0])[:, None]
    np.put_along_axis(ranks, order, np.broadcast_to(rows, order.shape), axis=0)
    u_a = ranks[:n_a].sum(axis=0) - n_a * (n_a - 1) / 2.0
    return 2.0 * u_a / (n_a * n_b) - 1.0


def continuous_null(tree: PhyloTree, idx_a: np.ndarray, idx_b: np.ndarray,
                    n_sims: int, rng) -> np.ndarray:
    """Null distribution of |Cliff's delta| between two tip sets under BM.

    Depends only on the tree and which tips are in which group -- **not on any
    feature's values** -- because a rank statistic is invariant to the scale
    and offset that are the only things a feature contributes under Brownian
    motion. So one call serves every feature sharing the same usable-tip set,
    which is what makes a large ``n_sims`` affordable here.
    """
    values = simulate_bm(tree, n_sims, rng)
    return np.abs(_cliffs_delta_columns(values, idx_a, idx_b))


def phylo_p_from_null(observed: float, null: np.ndarray) -> float:
    """``(count(null >= |observed|) + 1) / (n + 1)``, matching the permutation
    p-value convention in :func:`mpph.analysis.permanova`.

    The ``- 1e-12`` slack is load-bearing, not cosmetic: the observed
    statistic is a difference of small-integer ratios, and an exact ``>=``
    drops legitimate ties at the 1e-17 level -- on the clade-confounded
    acceptance case that is the difference between p ~ 0.8 and p ~ 0.0.
    """
    if null.size == 0:
        return float("nan")
    hits = int(np.sum(null >= abs(observed) - 1e-12))
    return (hits + 1) / (null.size + 1)


def binary_phylo_test(tree: PhyloTree, states: np.ndarray, idx_a: np.ndarray,
                      idx_b: np.ndarray, *, n_sims: int, rng,
                      min_accepted: int = 100) -> dict:
    """Phylogenetically-corrected p-value for a presence/absence feature.

    Simulates the feature on the tree under an ML-fitted symmetric Mk rate
    (a null with no group association by construction) and compares the
    observed prevalence difference to that null.

    **The null is conditioned on the observed number of present tips.** This
    is not a refinement, it is what makes the test work: without it roughly a
    third of simulated replicates come out invariant, contribute a zero
    statistic, dilute the tail, and the test then calls the pure
    clade-confounded case significant (measured p = 0.013 at 32 tips, getting
    worse as the tree grows -- i.e. it fails on exactly the scenario it exists
    to catch). Conditioning gives p ~ 0.8 there. Fisher's exact test
    conditions on its margins; this null has to as well.

    Acceptance uses exact ``k`` where it can, widening the band by one only
    while too few replicates are accepted and never past 10% of the tip
    count; the band actually used is returned so a reader can audit it.
    """
    used = np.concatenate([idx_a, idx_b])
    obs = states[used]
    k_obs = int(np.nansum(obs > 0.5))
    prev_a = float(np.nanmean(states[idx_a] > 0.5))
    prev_b = float(np.nanmean(states[idx_b] > 0.5))
    observed = prev_a - prev_b

    q_hat, at_boundary = fit_mk_rate(tree, states)
    sims = simulate_mk(tree, q_hat, n_sims, rng)
    k_sim = sims[used].sum(axis=0)

    tol_cap = max(0, int(np.floor(0.10 * len(used))))
    tol = 0
    while True:
        accept = np.abs(k_sim - k_obs) <= tol
        if int(accept.sum()) >= min_accepted or tol >= tol_cap:
            break
        tol += 1

    n_accepted = int(accept.sum())
    warning = ""
    if at_boundary:
        warning = "mk_rate_at_boundary"
    if n_accepted == 0:
        return {"p_value_phylo": float("nan"), "n_phylo_replicates_accepted": 0,
                "phylo_k_tolerance": tol, "mk_rate": q_hat,
                "phylo_warning": warning or "no_replicates_matched_prevalence"}

    kept = sims[:, accept]
    sim_prev_a = kept[idx_a].mean(axis=0)
    sim_prev_b = kept[idx_b].mean(axis=0)
    null = np.abs(sim_prev_a - sim_prev_b)
    if tol > 0 and not warning:
        warning = f"prevalence_band_widened_to_{tol}"
    return {"p_value_phylo": phylo_p_from_null(observed, null),
            "n_phylo_replicates_accepted": n_accepted,
            "phylo_k_tolerance": tol, "mk_rate": q_hat,
            "phylo_warning": warning}


def fitch_changes(tree: PhyloTree, states: np.ndarray) -> int:
    """Minimum number of 0/1 state changes on the tree (Fitch parsimony).

    Deterministic, O(nodes), no RNG, and exactly hand-checkable -- which is
    why it is the diagnostic shipped alongside the simulated p-value rather
    than a phylogenetic-signal index like Fritz & Purvis' D (noisy per
    feature, and its observed statistic ignores the branch lengths its own
    null uses).

    ``1`` means the trait arose exactly once. If that single origin sits on
    the branch that also separates the two groups being compared, then no
    method can distinguish association from coincidence -- the effective
    sample size for the comparison is 1 (Maddison & FitzJohn 2015). A large
    phylogenetic p-value there is the correct answer, not a broken test.

    ``states`` is indexed by tip order (``tree.tip_labels``); NaN means
    unknown and contributes no constraint.
    """
    sets: list = [None] * tree.n_nodes
    tips = tree.tip_indices
    for pos, node in enumerate(tips):
        value = states[pos]
        if np.isnan(value):
            sets[node] = {0, 1}
        else:
            sets[node] = {int(value > 0.5)}

    changes = 0
    for i in range(tree.n_nodes - 1, -1, -1):
        if tree.is_tip[i]:
            continue
        child_sets = [sets[c] for c in tree.children[i]]
        common = set.intersection(*child_sets) if child_sets else set()
        if common:
            sets[i] = common
        else:
            sets[i] = set.union(*child_sets)
            changes += 1
    return int(changes)
