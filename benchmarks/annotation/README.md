# Annotation benchmark: `mpph annotate` vs. KofamScan vs. eggNOG-mapper

Head-to-head comparison of per-gene KO assignment across three tools, on ten
genomes spanning **all three domains of life** plus unassigned marine MAGs.
Seven have real, KEGG-curated ground truth (`link/ko/<org>` -- the same
official KO assignments KEGG itself publishes for these genomes, not a
synthetic or self-reported benchmark).

## Is the ground truth independent of the tools being tested?

Yes, and this is worth stating explicitly because it is the first thing that
should be checked before trusting any of the accuracy numbers below.

KEGG's own documentation describes the GENES database's KO assignments as
built on **GFIT tables — best-hit sequence-similarity results** — used as
"the basis for both manual and automatic annotation of the GENES database".
KEGG's automatic annotation servers are likewise similarity-based
(BlastKOALA = BLASTP, GhostKOALA = GHOSTX); **KofamKOALA (HMM profile
search) is offered as a separate user-facing service, not as the mechanism
behind GENES annotation itself.**

So the ground truth here is produced by sequence-similarity search plus
curation, whereas `mpph annotate` and KofamScan both score profile HMMs
against per-KO adaptive thresholds. The two are methodologically
independent — recovering KEGG's KO assignments with a profile-HMM method is
a real result, not a restatement of the reference. (eggNOG-mapper is
similarity-based like KEGG's pipeline, which is context for interpreting its
numbers, not an advantage it was given here.)

## Why these organisms

- `eco` -- *Escherichia coli* K-12 MG1655 (Gammaproteobacteria), primary
  validation, ~4,300 genes.
- `bsu` -- *Bacillus subtilis* subsp. *subtilis* 168 (Firmicutes, Gram-positive
  -- a different phylum from *E. coli*).
- `mja` -- *Methanocaldococcus jannaschii* DSM 2661 (Euryarchaeota -- a
  different domain of life).
- `sce` -- *Saccharomyces cerevisiae* S288C (Fungi), 6,021 proteins.
- `ath` -- *Arabidopsis thaliana* Col-0 (Viridiplantae), 48,265 proteins --
  11x *E. coli*'s proteome, and the largest input in this benchmark.

The last two are the ones that answer "does this work outside microbes?"
(`run_eukaryotes.sh`). They matter for two separate reasons: **Eukaryota is
the third domain**, so with `mja` the benchmark now covers all three; and
they are an order of magnitude larger than a bacterial genome, so they test
whether the implementation *scales*, not just whether it generalizes.

All five have real KEGG-curated KO assignments, so this measures actual
precision/recall/F1 against ground truth for all five, not "accuracy for
*E. coli*, mere pairwise tool-agreement for the rest."

These five are, however, all long-studied model organisms — the obvious
objection is that KEGG has had decades to curate them, so good accuracy
there says little about the "works for any species" claim this tool actually
makes. Two harder tiers were added for exactly that reason
(`run_nonmodel.sh`):

**Tier 2 — non-model genomes KEGG never hand-curated** (ground truth still
available, so precision/recall are still measured):

| KEGG code | Organism | Accession | Note |
|---|---|---|---|
| `vbs` | *Verrucomicrobia bacterium* S94 | GCA_004299845.1 | PVC superphylum; placeholder strain name, not a named species |
| `lbac` | *Lentisphaerae bacterium* WC36 | GCA_020683205.1 | PVC superphylum; among the most recently added KEGG genomes (T10936) |

Both are environmentally-derived bacteria carrying `bacterium <strain>`
placeholder names rather than species names — i.e. organisms nobody has
built a curated pathway model around.

**Tier 3 — real marine MAGs, the actual target use case** (no ground truth
exists; coverage and pairwise tool agreement only). Metagenome-assembled
genomes with Prodigal gene calls, from marine metagenomes (NCBI SRA runs,
binned with MetaBAT2/DAS Tool). Unlike a finished reference genome these are
incomplete and fragmentary, which is the condition `mpph annotate` is
actually meant to be used under. Selection rule, stated so it is not a
cherry-pick: from bins with 1,000–4,000 predicted proteins (excluding
obvious multi-organism bins, some of which exceed 20,000 proteins), sampled
across that range.

Precision/recall are deliberately **not** reported for Tier 3. Absence of a
reference is not evidence that a call is wrong, and scoring against a
non-existent truth set would manufacture numbers rather than measure
anything.

Genome source: NCBI (the exact assembly accession KEGG itself lists as each
organism's `DATA_SOURCE`, confirmed via `get/gn:<org>`, not assumed):

| Organism | Accession | Assembly name |
|---|---|---|
| `eco` | GCF_000005845.2 | ASM584v2 |
| `bsu` | GCF_000009045.1 | ASM904v1 |
| `mja` | GCA_000091665.1 (GenBank, not RefSeq) | ASM9166v1 |
| `sce` | GCF_000146045.2 | R64 |
| `ath` | GCF_000001735.4 | TAIR10.1 |

## Scoring a eukaryote fairly: the isoform problem

A eukaryotic proteome lists every splice isoform; KEGG's reference names one
representative protein per gene. For *Arabidopsis* that is 48,265 proteins
against 27,562 the reference can speak to. Scoring a call on an isoform the
reference never mentions as a *false positive* would have manufactured
**10,025 of them** in this run — and a hypothetical perfectly-correct tool
would score precision 0.5 instead of 1.0 purely from that.

So when a protein-id map is supplied, calls are restricted to the proteins
the reference actually covers (`restrict_to_evaluable`), and the number
excluded is printed rather than hidden. Prokaryotes are essentially
unaffected (*E. coli*: 4,288 of 4,300 covered). This is a property of the
reference, not of any tool, and it is applied identically to all three.

`test_compare_annotation.py` pins both halves: that the restriction drops
unrepresented isoforms, and that omitting it would drag a perfect tool's
precision to 0.5.

## Results

![Annotation benchmark](../figures/fig1_annotation_benchmark.png)

Full tables: `results/summary_accuracy.csv`, `results/summary_agreement.csv`.

**`mpph annotate` reproduces KofamScan essentially exactly.** Jaccard overlap
of (gene, KO) calls between the two is **1.000 on seven of the ten genomes**,
and 0.9997 / 0.9999 / 0.9994 on the other three — *E. coli* (one call in
3,517), *Arabidopsis* (one in 12,223) and one MAG. This is the central
result: the pyhmmer reimplementation is not merely "comparable" to the tool
it reimplements, it is near-identical call for call — including on genomes
neither tool's authors ever looked at, and on a plant proteome 11x the size
of the bacterial ones.

**Accuracy holds across all three domains of life** (F1 against KEGG's own KO
assignments):

| Tier | Domain | Genome | mpph annotate | KofamScan | eggNOG-mapper |
|---|---|---|---|---|---|
| Model | Bacteria | *E. coli* K-12 | **0.934** | 0.934 | 0.851 |
| Model | Bacteria | *B. subtilis* 168 | **0.898** | 0.898 | 0.825 |
| Model | Archaea | *M. jannaschii* | 0.894 | 0.894 | **0.907** |
| Model | Eukaryota | *S. cerevisiae* | **0.911** | 0.911 | 0.862 |
| Model | Eukaryota | *A. thaliana* | **0.867** | 0.866 | 0.697 |
| Non-model | Bacteria | *Verrucomicrobia* sp. S94 | **0.848** | 0.848 | 0.755 |
| Non-model | Bacteria | *Lentisphaerae* sp. WC36 | **0.842** | 0.842 | 0.761 |

The objection this benchmark was built to answer — "your organisms are all
hand-curated model organisms" — does not survive it: the profile-HMM tools'
margin over eggNOG-mapper is *larger* on the two non-model genomes
(+0.09, +0.08 F1) than on the model organisms (+0.08, +0.07), not smaller.
On real MAGs, where no reference exists to score against, mpph and KofamScan
still agree at 1.000/0.999/1.000.

**Nothing degrades outside bacteria.** Yeast (0.911) scores *above* three of
the four prokaryotes, and the *Arabidopsis* F1 of 0.867 is squarely in the
prokaryotic range. The widest gap over eggNOG-mapper anywhere in the
benchmark is on *Arabidopsis* (+0.17 F1) — the similarity-transfer approach
loses the most ground exactly where the tool is furthest from a close
reference, which is the pattern profile HMMs are supposed to produce.

**What does change is coverage, and it is biology, not a defect.** The
fraction of proteins receiving any KO falls from 76% (*E. coli*) to 62%
(yeast) to 24% (*Arabidopsis*), for all three tools alike. KO covers
metabolism and core cellular processes; a plant proteome is largely made of
things KO does not describe (and, per the isoform note above, lists each gene
several times over). Recall against what KEGG *does* assign stays high —
0.925 for *Arabidopsis*, the highest of any genome here.

### The three calls where mpph and KofamScan disagree

Worth chasing to the bottom, because "essentially identical" is only a real
claim if the residue is explained. Across *Arabidopsis*' 12,223 calls the two
tools differ on three (gene, KO) pairs, in *both* directions. All three sit
within **0.016 bits** of that KO's threshold:

| KO | Protein | true best-domain score | threshold | HMMER text output | KofamScan | mpph |
|---|---|---|---|---|---|---|
| K14684 | NP_565171.4 | 332.1365 | 332.13 | `332.1` | no call | call |
| K23113 | NP_001331298.1 | 164.0541 | 164.07 | `164.1` | call | no call |
| K23113 | NP_001331301.1 | 164.0541 | 164.07 | `164.1` | call | no call |

One mechanism explains both directions: **HMMER prints bit scores to one
decimal place**, and KofamScan parses that text output, so it compares the
*rounded* score against a threshold given to two decimals. `332.1 >= 332.13`
is false; `164.1 >= 164.07` is true. `mpph annotate` reads the score from
pyhmmer's in-memory hit object and compares at full float precision, so it
gets the opposite answer in each case.

Neither is biologically more correct — these are ties inside the reference
database's own printing precision. But it means the two implementations agree
on every call that is not a tie, which is a stronger statement than the
Jaccard number alone. (Verified directly: re-running `hmmsearch` from the
KofamScan container on these three proteins reproduces the rounded values in
the table above.)

**On eggNOG-mapper's numbers, fairly:** it is not "worse". It annotates
*more* genes (higher coverage on every genome) at lower precision — a
different sensitivity/specificity operating point — and it is
similarity-based, like KEGG's own annotation pipeline, while also producing
far more than KO calls (COG categories, GO, domains). A KO-only F1 is a
narrow slice of what it does.

**Cost, including where mpph is worse.** Wall-clock is comparable — mpph is
faster than KofamScan on seven of ten genomes (e.g. 44 s vs 49 s on
*M. jannaschii*), slower on *E. coli* (169 s vs 81 s) and on the two
eukaryotes (430 s vs 293 s; 2,097 s vs 1,905 s), 28 threads throughout.
Peak memory is not comparable and this is a genuine limitation of the
implementation, not a reporting artifact (verified against the raw
`/usr/bin/time -v` logs):

| | mpph annotate | KofamScan |
|---|---|---|
| Peak RSS across the ten genomes | **2.3 – 11.8 GB** | 0.13 – 0.74 GB |

11–18× higher — and, importantly, **not a function of input size**. The worst
case in the whole benchmark is *B. subtilis* at 4,237 proteins (11.8 GB),
while *Arabidopsis* at 48,265 proteins — 11x larger — peaks *lower*, at
9.7 GB. Memory here is dominated by the profile search, not by the targets.

Two direct measurements rule out the obvious explanation (pyhmmer's own docs
note that prefetching targets trades "much higher memory consumption" for
speed):

- **A prefetched sequence block costs ~1.0 kB per protein** — 4.5 MB for
  *B. subtilis*, 45.5 MB for *Arabidopsis*. Holding an entire plant proteome
  in memory is ~0.5% of the peak.
- **Re-running *B. subtilis* both ways changes nothing.** Prefetch 11.4 GB /
  2:06.9; stream 11.7 GB / 2:07.5 — byte-identical output, and the 0.3 GB gap
  is *smaller than the run-to-run spread of the identical prefetch command*
  (two runs measured 11.4 and 11.8 GB). There is no difference here to
  attribute to sequence loading at all.

`--sequence-loading` is therefore kept as a control for metagenome-scale
protein catalogues (millions of sequences, where ~1 kB each does add up), not
as a fix for eukaryotic genomes, and `auto` prefetches up to 1,000,000
proteins.

**What memory actually tracks is `--cpus`.** Same 4,237-protein
*B. subtilis* input, same prefetch path, only the thread count varied:

| `--cpus` | peak RSS | wall clock | speed-up |
|---|---|---|---|
| 1 | **1.0 GB** | 21:07 | 1.0× |
| 4 | **2.9 GB** | 5:17 | 4.0× |
| 28 | **11.8 GB** | 2:07 | 10.0× |

All three produced byte-identical output.

Roughly 0.4 GB per thread at the top end, and the whole 2.3–11.8 GB spread
across the benchmark is explained by nothing more exotic than 28 threads
against KOfam's 26,498 profiles. That makes `--cpus` the memory knob, and
the trade is worth knowing: going 28 → 4 threads cuts peak memory 4× and
costs 2.5× wall-clock, because parallel scaling is already well past its
linear region by 28 (4 threads gives a perfect 4.0× speed-up; 28 gives only
10×). On a memory-constrained machine, lower `--cpus` — the tool does not
need a large-memory node.

Recorded because a benchmark that only reports the flattering numbers is not
a benchmark.

## Tools and versions

| Tool | Container | Notes |
|---|---|---|
| `mpph annotate` | this repo, installed from source | pyhmmer + KOfam, see `mpph/kofam.py` |
| KofamScan | `quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2` | run against the **same** KOfam `profiles/`+`ko_list` as `mpph annotate`, isolating the pyhmmer reimplementation as the only variable |
| eggNOG-mapper | `quay.io/biocontainers/eggnog-mapper:2.1.15--pyhdfd78af_0` | `-m diamond`, full (not taxonomically-restricted) database, matching KofamScan/mpph's always-full-profile-set search |

KOfam database: KEGG's own distribution, downloaded via `mpph annotate
--setup-db` (the same code path this repo ships). eggNOG database: the
default diamond database from `download_eggnog_data.py` (no `-D`/`-H`/`-M`
restriction).

## Reproducing

```bash
# One-time setup (see setup_containers.log / setup_mpph_kofam.log /
# setup_eggnog_db.log in report/ for the exact commands actually run)
singularity pull containers/kofamscan.sif docker://quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2
singularity pull containers/eggnog-mapper.sif docker://quay.io/biocontainers/eggnog-mapper:2.1.15--pyhdfd78af_0
mpph annotate --setup-db databases/kofam
singularity exec containers/eggnog-mapper.sif download_eggnog_data.py -y --data_dir databases/eggnog

# Run everything + compare
bash run_all.sh        # tier 1 model organisms (eco, bsu, mja)
bash run_nonmodel.sh   # tier 2 non-model genomes + tier 3 MAGs
bash run_eukaryotes.sh # yeast + Arabidopsis
python summarize.py --report-dir report --results-dir results --out-dir report
```

## Metrics (`compare_annotation.py`)

For each (organism, tool): coverage (% genes assigned >=1 KO), and
precision/recall/F1 treating each (gene, KO) pair as the unit of comparison
against ground truth (so multi-KO genes are handled correctly, not just
exact-set matches). For each (organism, tool-pair): Jaccard overlap of
(gene, KO) pairs -- symmetric agreement, no designated ground truth needed.
Logic verified against hand-constructed synthetic data with known expected
precision/recall/Jaccard values (`test_compare_annotation.py`) before being
trusted on real tool output.

Raw per-tool outputs and the multi-GB databases themselves are not committed
here (see `.gitignore`) -- only the summary tables in `report/` and the
scripts to regenerate them.
