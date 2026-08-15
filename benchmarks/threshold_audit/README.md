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

## Reproducing

```bash
bash run_audit.sh
python summarize_audit.py --audit-dir threshold_audit --ground-truth-dir ground_truth \
    --organisms eco bsu mja vbs lbac sce ath \
    --out-summary results/summary.csv --out-curve results/curve.csv
python calibration_transfer.py --audit-dir threshold_audit --ground-truth-dir ground_truth \
    --train eco bsu --test mja vbs lbac \
    --out results/transfer.csv --out-reliability results/reliability.csv
```

## Caveat that must travel with these numbers

KEGG's reference is itself similarity-derived and incomplete, so a hit counted
here as a false positive may be a genuine function KEGG has not recorded.
Measured precision is therefore a **lower bound**, not a point estimate. That
is the conservative direction for the argument, but the numbers must not be
requoted later as point estimates.
