"""Three loose ends around the arbitration result, none of them flattering.

**1. The eggNOG-mapper comparison is currently unfair, in this tool's favour.**
eggNOG-mapper also emits several KOs for one gene (15.4% of *E. coli* genes,
9.9% of *Arabidopsis*), so comparing an arbitrated `mpph annotate` against an
unarbitrated eggNOG measures the arbitration step, not the two annotators. It
cannot be arbitrated the *same* way -- its `KEGG_ko` column carries no per-KO
score to rank by -- which is a genuine structural advantage of the profile-HMM
route and should be stated as one rather than silently enjoyed. Reported here:
eggNOG scored (a) as published, and (b) under a first-listed-KO-per-gene rule,
which is the best like-for-like bound available without scores.

**2. The per-genome gains were published as point estimates.** Paired
bootstrap CIs over calls, plus the across-genome sign test.

**3. Which KOs lose arbitration most often**, and do they concentrate in
particular families -- the mechanism check for the recall-loss bias result.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "annotation"))
from compare_annotation import load_ground_truth, load_protein_id_map

_KO = re.compile(r"K\d{5}")


def _strip(gene: str) -> str:
    return gene.rsplit(".", 1)[0] if gene.rsplit(".", 1)[-1].isdigit() else gene


def emapper_pairs(path: Path, *, first_only: bool) -> set[tuple[str, str]]:
    """(gene, KO) pairs from an emapper .annotations file.

    ``first_only`` keeps just the first KO listed per gene. eggNOG gives no
    per-KO score, so "first listed" is the only one-per-gene rule available --
    it is a *bound*, not an equivalent of margin-based arbitration, and is
    labelled as such wherever it is reported.
    """
    header = None
    idx_ko = None
    out: set[tuple[str, str]] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("##") or not line.strip():
            continue
        f = line.lstrip("#").split("\t")
        if header is None:
            header = f
            idx_ko = header.index("KEGG_ko")
            continue
        if line.startswith("#"):
            continue
        kos = _KO.findall(f[idx_ko]) if idx_ko < len(f) else []
        if not kos:
            continue
        gene = _strip(f[0])
        out.update((gene, k) for k in (kos[:1] if first_only else kos))
    return out


def mapper_pairs(path: Path) -> set[tuple[str, str]]:
    out = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        f = line.split("\t")
        if len(f) >= 2 and f[1].startswith("K"):
            out.add((_strip(f[0]), f[1]))
    return out


def prf(pred: set, truth: set) -> dict:
    tp = len(pred & truth)
    precision = tp / len(pred) if pred else float("nan")
    recall = tp / len(truth) if truth else float("nan")
    f1 = 2 * tp / (len(pred) + len(truth)) if (pred or truth) else float("nan")
    return {"n_calls": len(pred), "precision": precision, "recall": recall,
            "f1": f1}


def paired_bootstrap(pred_a: set, pred_b: set, truth: set, *,
                     n_boot: int = 2000, seed: int = 0) -> dict:
    """CI for (b - a) in precision and F1, resampling *genes* not calls.

    Genes are the independent unit here -- calls on the same gene are exactly
    what arbitration acts on, so resampling calls would break the dependency
    the statistic is about and give CIs that are too narrow.
    """
    genes = sorted({g for g, _ in pred_a | pred_b | truth})
    rng = np.random.default_rng(seed)
    by_gene_a, by_gene_b, by_gene_t = {}, {}, {}
    for g, k in pred_a:
        by_gene_a.setdefault(g, set()).add((g, k))
    for g, k in pred_b:
        by_gene_b.setdefault(g, set()).add((g, k))
    for g, k in truth:
        by_gene_t.setdefault(g, set()).add((g, k))

    d_p, d_f = np.empty(n_boot), np.empty(n_boot)
    for i in range(n_boot):
        pick = rng.choice(len(genes), len(genes))
        a, b, t = set(), set(), set()
        for j, gi in enumerate(pick):
            g = genes[gi]
            tag = f"#{j}"          # keep resampled copies distinct
            a |= {(g + tag, k) for _, k in by_gene_a.get(g, ())}
            b |= {(g + tag, k) for _, k in by_gene_b.get(g, ())}
            t |= {(g + tag, k) for _, k in by_gene_t.get(g, ())}
        ra, rb = prf(a, t), prf(b, t)
        d_p[i] = rb["precision"] - ra["precision"]
        d_f[i] = rb["f1"] - ra["f1"]
    return {"d_precision_lo": float(np.percentile(d_p, 2.5)),
            "d_precision_hi": float(np.percentile(d_p, 97.5)),
            "d_f1_lo": float(np.percentile(d_f, 2.5)),
            "d_f1_hi": float(np.percentile(d_f, 97.5))}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--ground-truth-dir", required=True)
    p.add_argument("--organisms", nargs="+", required=True)
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--out-fairness", required=True)
    p.add_argument("--out-ci", required=True)
    p.add_argument("--out-kos", required=True)
    args = p.parse_args()

    results, gt = Path(args.results_dir), Path(args.ground_truth_dir)
    fair_rows, ci_rows = [], []
    lost_kos: Counter = Counter()
    won_kos: Counter = Counter()

    for org in args.organisms:
        id_map_path = gt / f"protein_id_map_{org}.tsv"
        id_map = (load_protein_id_map(id_map_path, org)
                  if id_map_path.exists() else None)
        truth_map = load_ground_truth(gt / f"link_ko_{org}.tsv", org, id_map)
        truth = {(g, k) for g, kos in truth_map.items() for k in kos}
        evaluable = set(id_map.values()) if id_map else None

        # Bound as a default argument, not captured: `evaluable` is a loop
        # variable, and a late-binding closure here would silently apply the
        # wrong genome's evaluable set the moment this is called lazily.
        def restrict(s, evaluable=evaluable):
            return {(g, k) for g, k in s if g in evaluable} if evaluable else s

        all_p = restrict(mapper_pairs(results / f"mpph_{org}.tsv"))
        best_p = restrict(mapper_pairs(results / f"{org}_best.tsv"))
        emap = results / f"emapper_{org}.emapper.annotations"

        row = {"organism": org}
        for label, pred in (("mpph_all", all_p), ("mpph_best", best_p)):
            fair_rows.append({"organism": org, "tool": label, **prf(pred, truth)})
        if emap.exists():
            for label, first in (("eggnog_published", False),
                                 ("eggnog_first_ko_bound", True)):
                pred = restrict(emapper_pairs(emap, first_only=first))
                fair_rows.append({"organism": org, "tool": label,
                                  **prf(pred, truth)})

        ci = paired_bootstrap(all_p, best_p, truth, n_boot=args.n_boot)
        ci_rows.append({"organism": org,
                        "d_precision": prf(best_p, truth)["precision"]
                        - prf(all_p, truth)["precision"],
                        "d_f1": prf(best_p, truth)["f1"] - prf(all_p, truth)["f1"],
                        **ci})
        row.update(ci)

        lost_kos.update(k for _, k in (all_p & truth) - best_p)
        won_kos.update(k for _, k in best_p - all_p)

    fair = pd.DataFrame(fair_rows)
    ci = pd.DataFrame(ci_rows)
    kos = pd.DataFrame(lost_kos.most_common(30), columns=["ko", "times_lost"])

    for path, frame in ((args.out_fairness, fair), (args.out_ci, ci),
                        (args.out_kos, kos)):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, index=False)

    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print("=== 1. like-for-like comparison (eggNOG cannot be margin-arbitrated) ===")
        print(fair.round(4).to_string(index=False))
        print("\nmean by tool:")
        print(fair.groupby("tool")[["precision", "recall", "f1"]].mean()
              .round(4).to_string())
        print("\n=== 2. paired bootstrap CIs, resampling genes ===")
        print(ci.round(4).to_string(index=False))
        n_up = int((ci["d_precision"] > 0).sum())
        n = len(ci)
        # one-sided sign test, H0: improvement equally likely either way
        p_sign = 0.5 ** n
        print(f"\nprecision improved on {n_up}/{n} genomes; "
              f"one-sided sign test p = {p_sign:.4f}")
        print("\n=== 3. KOs most often losing arbitration ===")
        print(kos.head(15).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
