# Literature matrix — comparative metabolic inference from genomes

Compiled 2026-08-15. Purpose: establish, with citations, what already exists,
so that any claim of novelty for MPPH is bounded by evidence rather than by
not having looked. Every "does NOT" below is a checked property of the cited
work, not an assumption.

---

## A. Metabolic reconstruction / pathway profiling tools (the direct field)

| Tool | Venue | What it does | What it does **not** do |
|---|---|---|---|
| **DRAM** | Nucleic Acids Res 2020 | Annotates MAGs/viral contigs against KEGG, UniRef90, Pfam, dbCAN, MEROPS; "distills" into functional categories | No uncertainty on any call; no phylogenetic model; presence is a hard call |
| **anvi'o** `anvi-estimate-metabolism` | Genome Biol 2021 / ongoing | Pathwise + stepwise KEGG module completeness; `anvi-compute-metabolic-enrichment` for group comparison | **Hard 0.75 completion threshold**; no confidence interval on completeness; enrichment test treats genomes as independent |
| **MicrobeAnnotator** | BMC Bioinformatics 2020 | Iterative multi-DB annotation → KO → module completeness | Binary KO calls; no uncertainty; no phylogeny |
| **KEGG-Decoder** | — | KO list → pathway completeness heatmap | Fractions only; no statistics at all |
| **METABOLIC** | Microbiome 2022 | HMM-based metabolic/biogeochemical profiling, community scale | No uncertainty; no phylogenetic correction |
| **MetaPathPredict** | **eLife 2024** | **Deep learning predicts complete KEGG modules in *incomplete* genomes.** Trained on ~30,596 RefSeq+GTDB genomes, 305,960 obs incl. artificial down-sampling to 10–90% of genes. Mean **F1 0.96** at 30–90% completeness | **(1)** Outputs **binary** presence/absence, **not calibrated probabilities**. **(2)** Uses **binary KO vectors** — discards E-values/bit scores entirely. **(3)** **No phylogenetic control** on the train/test split, over a database (RefSeq) with extreme clade redundancy. Authors themselves report an "over-exuberant positive class prediction problem" below 30% |
| **MPPH** (this repo) | unpublished | KO annotation (≡ KofamScan), presence/completeness matrices, ORA/GSEA (≡ clusterProfiler), pathway maps, phylogeny-aware differential test | See §D — most components reproduce an existing tool by design |

**Field-level property:** every tool above converts evidence into a hard
call — a KO is present or absent, a module is complete or not — and every
downstream comparison then treats those hard calls as observed data.

---

## B. Phylogeny-aware association testing (the methods MPPH's `--tree` overlaps)

| Method | Venue | Approach |
|---|---|---|
| **Scoary** | Genome Biol 2016 | Pan-GWAS on accessory-gene presence/absence; pairwise-comparisons method counts minimum independent co-emergences on the tree |
| **treeWAS** | PLOS Comput Biol 2018 | Ancestral state reconstruction + simulation under a homoplasy distribution; three association tests |
| **pyseer** | Bioinformatics 2018 | k-mer/unitig/gene GWAS with population-structure control (mixed models, lineage effects) |
| **hogwash** | Microb Genom 2020 | Two ancestral-reconstruction methods (phyC + a stricter variant) for binary/continuous phenotypes |
| **phyC** | — | Phylogenetic convergence: mutations arising independently more often than chance |

**This is the honest problem for MPPH.** Correcting a gene/feature
presence–absence association for phylogenetic non-independence is a solved,
well-cited idea with at least four mature implementations. MPPH's
`compare --tree` applies it to KEGG features rather than accessory genes,
which is a **port, not an invention**. It is good engineering and it fixes a
real defect in *this* toolkit; it is not by itself a Nature Methods claim,
and the manuscript must not imply otherwise.

---

## C. The two papers that define the open problem

### C1. Genome incompleteness confounds functional inference
*Impact of microbial genome completeness on metagenomic functional inference*,
ISME Communications 3:12 (2023).

- Completeness 70% → 100% is associated with a **15 ± 10% increase in module
  fullness**. Missing genome ⇒ systematically understated metabolism.
- **The bias persists in "high-quality" (>90% complete) MAGs** — the standard
  CheckM filter does not remove it.
- **Phylum-dependent**: strongest in Proteobacteria, then Firmicutes,
  Actinobacteriota, Bacteroidota. So the distortion is itself phylogenetically
  structured — it is *confounded with* clade membership.
- **Module-size dependent**: modules with fewer steps are hit hardest.
- Their own correction "tended to **overcorrect**"; they state there is
  "ample room for improvement".

### C2. Population structure inflates ML accuracy in bacterial genomics
Yu et al., *Biased sampling driven by bacterial population structure confounds
machine learning prediction of antimicrobial resistance*, PLOS Biology (2025).

- >24,000 genomes, 5 species, 27 antibiotics, meta-analysis of 6,740 models.
- Models trained on realistic (biased) samples **conflate phylogenetic markers
  with genuine causal determinants**; AUC drops substantially under clade-held-out
  evaluation vs random splits.
- Adding more clades did **not** fix it — bias persisted.
- Explicit recommendation: **phylogeny-aware cross-validation, testing on
  held-out clades**, is required to demonstrate generalisability.

**The join nobody has made:** C2's critique has been made for AMR. It has
**not** been applied to metabolic-module prediction — where MetaPathPredict
reports F1 0.96 from random splits over a phylogenetically redundant database,
exactly the design C2 shows is inflated.

---

## D. What is genuinely unoccupied

Three distinct causes produce the *same observable* — a missing KO — and no
tool in §A distinguishes them:

1. **Detection**: the HMM bit score fell below KOfam's adaptive threshold. A
   score of 331.9 against a threshold of 332.13 is not the same evidence as no
   hit at all — but both are written as `0`. Every tool discards the score.
2. **Recovery**: the genome is incomplete. Per C1 this is large, persists past
   the 90% filter, and is phylum-structured.
3. **Evolution**: the capability is genuinely absent in this lineage.

Unoccupied, in order of how defensible the novelty claim is:

- **D1. Calibration.** No tool in §A reports a reliability diagram for its
  presence calls. MetaPathPredict reports F1, not calibration. "How often is a
  call labelled 90%-confident actually right?" is unanswered field-wide.
- **D2. Clade-held-out re-evaluation of metabolic prediction.** Directly
  testable, directly follows C2, directly targets the one recent high-profile
  competitor. A negative result (accuracy holds up) is publishable too.
- **D3. Joint separation of the three channels above** into one posterior,
  with the phylogenetic prior doing the work of distinguishing "assembly gap"
  from "true loss" — which is precisely the information §B's methods use, but
  applied to *imputation* rather than to *association*.

---

## E. Bar check — what Nature Methods actually publishes here

NM wants a method enabling something previously not possible, benchmarked
against the field's incumbents, with the benchmark itself credible. It does
not publish integration/convenience toolkits, however well engineered.

- Reproducing KofamScan call-for-call: **excellent evidence of correctness,
  zero novelty**. Belongs in validation, not in the claim.
- Matching clusterProfiler to 3 decimals: same.
- `compare --tree`: fixes a real defect, but see §B — prior art.
- **A calibrated, evidence-resolved posterior over metabolic capability,
  validated by clade-held-out evaluation and reliability diagrams, that
  demonstrably changes published conclusions**: this is the only framing in
  reach that is a *methods* contribution rather than a *tool* contribution.

Realistic venue ladder if the pilot in §D2 succeeds: Nature Methods is a
stretch but not absurd; Genome Biology / Nature Communications / Microbiome
are the likely landing zones. That should be said plainly rather than
promised away.
