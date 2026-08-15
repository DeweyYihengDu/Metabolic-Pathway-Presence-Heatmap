# Threshold audit: what the binary KO call throws away

Every tool in this space — KofamScan, `mpph annotate`, DRAM, MicrobeAnnotator,
anvi'o — turns the HMM search into a binary call at each KO's adaptive
threshold and discards the score. This measures what that costs, against
KEGG's own curated `link/ko/<org>` assignments.

Full write-up and interpretation: [`paper/nmeth/pilot_result_threshold_audit.md`](../../paper/nmeth/pilot_result_threshold_audit.md).

## Scripts

| Script | Purpose |
|---|---|
| `collect_subthreshold.py` | Re-runs the KOfam search retaining hits *below* threshold, emitting `delta = score − threshold` per (gene, KO) |
| `score_precision.py` | Precision vs `delta`, with a label-shuffled control for the background rate |
| `summarize_audit.py` | Cross-genome summary: the accepted-call calibration curve and the share of errors contributed by low-margin calls |
| `calibration_transfer.py` | **The gate.** Fits P(correct \| margin) on one set of genomes, tests calibration on a held-out lineage against a constant baseline, with paired bootstrap CIs |
| `run_audit.sh` | Driver across every genome with KEGG ground truth |

Logic is verified against hand-constructed data with known expected values
(`test_score_precision.py`, 23 tests) before being trusted on real output —
including that the truth join is *pairwise* rather than per-gene, which if
wrong would inflate the headline number.

## Results

`results/summary.csv`, `results/curve.csv`, `results/transfer.csv`,
`results/reliability.csv`.

**A KO call's precision spans 0.27 → 0.99 as a function of its margin above
threshold.** Every tool writes all of them as the same `1`.

| genome | domain | precision, all calls | % calls within 20 bits | precision of those | % of all false positives they carry |
|---|---|---|---|---|---|
| *E. coli* K-12 | Bacteria | 0.900 | 5.7% | 0.572 | 24.4% |
| *B. subtilis* 168 | Bacteria | 0.858 | 7.9% | 0.449 | 30.8% |
| *M. jannaschii* | Archaea | 0.903 | 12.9% | 0.652 | 46.5% |
| *Verrucomicrobia* sp. S94 | Bacteria (PVC) | 0.845 | 12.7% | 0.431 | 46.9% |
| *Lentisphaerae* sp. WC36 | Bacteria (PVC) | 0.848 | 11.3% | 0.597 | 29.9% |
| *S. cerevisiae* | Eukaryota | 0.945 | 5.6% | 0.752 | 25.4% |
| *A. thaliana* | Eukaryota | 0.815 | 7.6% | 0.319 | 28.0% |

6–13% of calls carry 24–47% of every error — a **2.7–4.5× enrichment in every
one of seven genomes spanning all three domains**. The low-margin fraction
roughly doubles from the model bacteria to the archaeon and the two
environmental PVC genomes, i.e. it is worst on exactly the genomes these tools
are built for. *A. thaliana* is the sharpest case: its low-margin calls are
correct less than a third of the time.

**The sub-threshold band is not a reservoir of lost signal**: precision below
the cut is ~0.10 in the first 5 bits and reaches background by −20. KOfam's
threshold is well placed. The information is in the *accepted* calls.

**Transfer (the gate):** fitted on *E. coli* + *B. subtilis*, tested on held-out
lineages, the margin model beats a constant-confidence baseline on Brier
significantly for both PVC genomes (+0.030, +0.015) — but never improves ECE,
and is significantly worse on the archaeon. **The ranking transfers; the
absolute probabilities do not.** A single global confidence curve is therefore
not shippable, and calibration has to be lineage-aware. That is a result about
what the method must do, not a failed experiment.

## The fix that came out of it: KO arbitration

The audit said the errors live on genes that received more than one KO. That
turns out to be a structural gap, not noise. **KOfam sets each KO's threshold
independently**, by maximising that KO's own F-measure in isolation — nothing
in the procedure makes KOs compete. A protein matching several related
profiles clears all of their thresholds at once, and KofamScan (and `mpph
annotate` by default, since it reproduces KofamScan) emits all of them.

KEGG's reference does not work that way. Across all seven benchmark genomes,
**≥99.88% of genes carry exactly one KO** (max observed: 2). So the extras are
over-calls by construction — on *E. coli* they are 13.5% of calls but **69% of
all false positives**.

`--multi-ko-policy best` adds the missing arbitration step: per gene, keep the
KO furthest above *its own* threshold. Ranking by margin rather than raw score
matters, because thresholds span ~30 to >2000 bits and raw scores are not
comparable between KOs.

| genome | precision (all → best) | F1 (all → best) |
|---|---|---|
| *E. coli* K-12 | 0.900 → **0.937** | 0.935 → **0.941** |
| *B. subtilis* 168 | 0.858 → **0.908** | 0.898 → **0.911** |
| *M. jannaschii* | 0.903 → **0.934** | 0.894 → **0.904** |
| *Verrucomicrobia* sp. S94 | 0.845 → **0.878** | 0.848 → **0.861** |
| *Lentisphaerae* sp. WC36 | 0.848 → **0.864** | 0.842 → **0.849** |
| *S. cerevisiae* | 0.945 → **0.960** | 0.911 → **0.914** |
| *A. thaliana* | 0.815 → **0.837** | 0.867 → **0.874** |
| **mean** | **0.874 → 0.903** | **0.885 → 0.893** |

Precision rises on all seven and F1 rises on all seven, so this is a net gain
rather than recall traded for precision. On *E. coli* it puts `mpph annotate`
ahead of the tool it reimplements: F1 0.941 vs KofamScan's 0.935.

**Rule choice was tested, not assumed.** `results/arbitration.csv` compares
five rules. Ranking by *relative* margin (score/threshold) is consistently
worse than absolute margin, contradicting the obvious scale argument. A tie
window and a minimum-gap gate both leave F1 flat (0.8927–0.8934 across every
parameter value tried), so the gain comes from arbitrating at all, not from
tuning — which is the reassuring outcome. `--min-ko-gap` is exposed for
callers who want to push precision further (mean 0.909 at 40 bits) at a known
recall cost.

**The default remains `all`**, so `mpph annotate` still reproduces KofamScan
call-for-call out of the box and the equivalence result stands. `best` is
opt-in.

### Pushing precision further: `--min-margin`

Dropping calls close to their threshold trades recall for precision along a
measured frontier (mean over the seven genomes):

| `--min-margin` | precision, no arbitration | precision **+ arbitration** | recall (+arb) | F1 (+arb) |
|---|---|---|---|---|
| 0 | 0.874 | **0.903** | 0.886 | 0.893 |
| 10 | 0.894 | **0.915** | 0.864 | 0.888 |
| 20 | 0.907 | **0.924** | 0.840 | 0.879 |
| 50 | 0.932 | **0.943** | 0.739 | 0.826 |
| 80 | 0.947 | **0.956** | 0.630 | 0.755 |

Two things this shows, both worth stating:

1. **Arbitration is not "being stricter".** Arbitration alone (precision
   0.903, recall 0.886) beats a raised margin alone at 10 bits (0.894, 0.876)
   on *both* axes, and the two stack at every point on the frontier.
   Competing-KO error and weak-hit error are different failures.
2. **F1 falls monotonically with `--min-margin`.** It is only worth using when
   a false positive genuinely costs more than a false negative — a property of
   the downstream analysis, not of the tool. The default is 0.

Verified end-to-end: `--multi-ko-policy best --min-margin 20` on
*M. jannaschii* gives precision 0.955 / recall 0.796, against KofamScan's
0.903 / 0.884 on the same input.

**Caveat that constrains the claim.** KEGG's one-KO-per-gene structure may be
partly a curation convention rather than pure biology, so some of this gain is
agreement with how the reference is built. What is not convention: KOfam
emitting several mutually exclusive orthology assignments for one protein is a
real artefact of per-KO threshold fitting, and those extra calls are wrong far
more often than the calls they accompany (0.39 vs 0.96 precision on *E. coli*).

## Reproducing

```bash
bash run_audit.sh
python summarize_audit.py --audit-dir threshold_audit --ground-truth-dir ground_truth \
    --organisms eco bsu mja vbs lbac sce ath \
    --out-summary results/summary.csv --out-curve results/curve.csv
python calibration_transfer.py --audit-dir threshold_audit --ground-truth-dir ground_truth \
    --train eco bsu --test mja vbs lbac \
    --out results/transfer.csv --out-reliability results/reliability.csv
python arbitrate.py --audit-dir threshold_audit --ground-truth-dir ground_truth \
    --organisms eco bsu mja vbs lbac sce ath --out results/arbitration.csv
```

The shipped implementation was checked against this offline analysis
end-to-end on three genomes (`eco`, `mja`, `vbs`): running
`mpph annotate --multi-ko-policy best` and scoring the real output through
`compare_annotation.py` reproduces `arbitration.csv`'s `top_delta` row to six
decimal places on precision, recall and F1 in all three cases.

## Caveat that must travel with these numbers

KEGG's reference is itself similarity-derived and incomplete, so a hit counted
here as a false positive may be a genuine function KEGG has not recorded.
Measured precision is therefore a **lower bound**, not a point estimate. That
is the conservative direction for the argument, but the numbers must not be
requoted later as point estimates.
