"""Does genome-based module prediction survive phylogeny-aware evaluation?

MetaPathPredict (eLife 2024) predicts whether a KEGG module is complete in an
incomplete genome, from that genome's KO vector, and reports mean F1 0.96. It
was trained and tested on ~30,596 RefSeq+GTDB genomes split **at random**.
That database is dominated by many strains of the same species and many
species of the same genus, so a random split routinely puts near-identical
relatives on both sides. Yu et al. (PLOS Biology 2025) showed this inflates
accuracy badly for AMR prediction and recommended clade-held-out evaluation;
it has not been applied to metabolic-module prediction.

Protocol, following MetaPathPredict's own design:

1. Full KO set per organism -> KEGG module completeness -> label = complete.
2. Down-sample the KO set to simulate an incomplete genome -> features.
3. Train to predict the label from the down-sampled vector.
4. Evaluate twice: random split, and **genus-blocked** split where no genus
   appears on both sides.

**The mandatory control.** Blocking reduces effective training diversity, so
some of any drop is expected from that alone rather than from leakage. Every
blocked split is therefore compared against a random split *of the same
training-set size*. Without this the experiment cannot distinguish "leakage
inflated the number" from "smaller training set performs worse", and the whole
result would be uninterpretable.

The model is deliberately a plain gradient-boosted tree rather than a
reimplementation of MetaPathPredict's network: the claim being tested is about
the *evaluation protocol*, not about any particular architecture, and the
comparison is between split schemes on one fixed model.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import f1_score
from sklearn.model_selection import GroupShuffleSplit, ShuffleSplit

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from mpph.kegg import DEFAULT_CACHE, make_session
from mpph.modules import fetch_module_definitions, list_modules, module_completeness


def load_ko_sets(ko_dir: Path) -> dict[str, set[str]]:
    out = {}
    for f in sorted(ko_dir.glob("*.tsv")):
        kos = set()
        for line in f.read_text(encoding="utf-8").splitlines():
            parts = line.split("\t")
            if len(parts) >= 2:
                kos.add(parts[1].split(":", 1)[-1])
        if kos:
            out[f.stem] = kos
    return out


def select_balanced_modules(ko_sets, defs, n_modules, rng, sample=250):
    """Pick usably-balanced modules from a *sample* of organisms.

    Scoring every module against every organism is 500+ x 3000 evaluations of
    a recursive boolean definition, which dominates the whole run. Prevalence
    only needs to be estimated well enough to reject the near-constant
    modules, and a few hundred organisms does that -- so the expensive full
    labelling happens once, for the modules actually used.
    """
    orgs = sorted(ko_sets)
    probe = [orgs[i] for i in rng.choice(len(orgs),
                                         size=min(sample, len(orgs)),
                                         replace=False)]
    keep = []
    for n, (mid, (definition, _, _)) in enumerate(defs.items(), 1):
        vals = [module_completeness(definition, ko_sets[o], defs) >= 1.0
                for o in probe]
        prev = float(np.mean(vals))
        if 0.15 <= prev <= 0.85:
            keep.append((mid, abs(0.5 - prev)))
        if n % 100 == 0:
            print(f"    screened {n}/{len(defs)} modules, {len(keep)} balanced",
                  flush=True)
    keep.sort(key=lambda t: t[1])
    return [m for m, _ in keep[:n_modules]]


def build_dataset(ko_sets, defs, module_ids, rng, retain=(0.3, 0.9)):
    """Features: down-sampled KO vector. Labels: module complete in the FULL set.

    Down-sampling is per organism with a retention fraction drawn uniformly in
    `retain`, matching MetaPathPredict's 10-90% simulated incompleteness. The
    label comes from the *complete* KO set -- the task is recovering what the
    genome really encodes from a partial view of it.
    """
    orgs = sorted(ko_sets)
    all_kos = sorted({k for s in ko_sets.values() for k in s})
    ko_index = {k: i for i, k in enumerate(all_kos)}
    X = np.zeros((len(orgs), len(all_kos)), dtype=np.int8)
    labels = {m: np.zeros(len(orgs), dtype=np.int8) for m in module_ids}
    for i, org in enumerate(orgs):
        full = ko_sets[org]
        keep_n = int(len(full) * rng.uniform(*retain))
        kept = rng.choice(sorted(full), size=keep_n, replace=False)
        X[i, [ko_index[k] for k in kept]] = 1
        for mid in module_ids:
            if module_completeness(defs[mid][0], full, defs) >= 1.0:
                labels[mid][i] = 1
        if (i + 1) % 500 == 0:
            print(f"    labelled {i + 1}/{len(orgs)} organisms", flush=True)
    return orgs, X, labels


def evaluate(X, y, groups, *, blocked: bool, train_size: int, seed: int) -> float:
    """Mean F1 over 3 splits at a fixed training-set size."""
    splitter = (GroupShuffleSplit if blocked else ShuffleSplit)
    scores = []
    test_frac = 0.25
    sp = splitter(n_splits=3, test_size=test_frac, random_state=seed)
    args = (X, y, groups) if blocked else (X, y)
    for tr, te in sp.split(*args):
        if train_size < len(tr):
            rs = np.random.default_rng(seed)
            tr = rs.choice(tr, size=train_size, replace=False)
        if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
            continue
        model = HistGradientBoostingClassifier(max_iter=120, random_state=seed)
        model.fit(X[tr], y[tr])
        scores.append(f1_score(y[te], model.predict(X[te]), zero_division=0))
    return float(np.mean(scores)) if scores else float("nan")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir", required=True)
    p.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    p.add_argument("--n-modules", type=int, default=40)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    data = Path(args.data_dir)
    with open(data / "taxonomy.tsv", encoding="utf-8") as fh:
        tax = {r["org"]: r for r in csv.DictReader(fh, delimiter="\t")}
    ko_sets = load_ko_sets(data / "ko")
    ko_sets = {o: s for o, s in ko_sets.items() if o in tax and tax[o]["genus"]}
    print(f"{len(ko_sets)} organisms with KO sets and a genus")

    session = make_session()
    cache = Path(args.cache_dir)
    defs = fetch_module_definitions(session, sorted(list_modules(session, cache)),
                                    cache)
    defs = {m: d for m, d in defs.items() if d[2] == "Pathway"}
    print(f"{len(defs)} Pathway modules")

    rng = np.random.default_rng(args.seed)
    # Modules with a usable class balance -- one present in 99% of genomes is
    # predicted perfectly by a constant and says nothing about either split
    # scheme. Screened on a sample first; see select_balanced_modules.
    print("screening modules for class balance ...", flush=True)
    module_ids = select_balanced_modules(ko_sets, defs, args.n_modules, rng)
    print(f"{len(module_ids)} balanced modules selected")

    print("building feature matrix and labels ...", flush=True)
    orgs, X, labels = build_dataset(ko_sets, defs, module_ids, rng)
    groups = np.array([tax[o]["genus"] for o in orgs])
    print(f"feature matrix {X.shape}, {len(set(groups))} genera")

    rows = []
    for mid in module_ids:
        y = labels[mid]
        if len(np.unique(y)) < 2:
            continue
        # Size-matched: the blocked split's training size is the cap for both.
        gss = GroupShuffleSplit(n_splits=1, test_size=0.25,
                                random_state=args.seed)
        tr_idx, _ = next(gss.split(X, y, groups))
        n_train = len(tr_idx)
        f1_blocked = evaluate(X, y, groups, blocked=True,
                              train_size=n_train, seed=args.seed)
        f1_random = evaluate(X, y, groups, blocked=False,
                             train_size=n_train, seed=args.seed)
        rows.append({"module": mid, "prevalence": float(y.mean()),
                     "n_train": n_train,
                     "f1_random_split": f1_random,
                     "f1_genus_blocked": f1_blocked,
                     "inflation": f1_random - f1_blocked})
        print(f"  {mid} prev={y.mean():.2f} random={f1_random:.3f} "
              f"blocked={f1_blocked:.3f} inflation={f1_random - f1_blocked:+.3f}",
              flush=True)

    df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    d = df.dropna(subset=["inflation"])
    print(f"\n=== {len(d)} modules evaluated, training size matched ===")
    print(f"mean F1, random split   : {d['f1_random_split'].mean():.4f}")
    print(f"mean F1, genus-blocked  : {d['f1_genus_blocked'].mean():.4f}")
    print(f"mean inflation          : {d['inflation'].mean():+.4f}")
    print(f"modules inflated        : {(d['inflation'] > 0).sum()} of {len(d)}")
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
