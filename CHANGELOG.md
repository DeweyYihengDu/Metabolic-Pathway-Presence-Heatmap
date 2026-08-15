# Changelog

## 3.21.1

**Correction to how 3.20.0 and 3.21.0 were framed: the precision gain does not
reach the module level.** Those releases reported `--multi-ko-policy best` as
an improvement measured on per-gene (gene, KO) pairs, which it is. What was
never measured, and should have been before shipping the claim, is whether it
improves what this toolkit actually reports downstream. Measured now, it does
not:

| level of the analysis | precision | recall | F1 |
|---|---|---|---|
| per-gene (gene, KO) pairs | +2.9 pts | -1.3 | **+0.8** |
| KO set (deduplicated, what completeness consumes) | +1.6 pts | -1.4 | **+0.0006** |

Against module completeness computed from KEGG's own KO sets, over seven
genomes and 522 Pathway modules: mean absolute error unchanged (0.0268 ->
0.0266), signed bias 44% **more** negative (-0.0122 -> -0.0176), and modules
complete in truth but broken in the estimate **up from 13.7 to 15.7** per
genome.

The recall arbitration gives up is also not uniform. Of 333 true calls
dropped, 17 pathways survive BH correction against a 2.78% baseline -- but
eight of those are one signal counted repeatedly, since KEGG's
neurodegeneration and metabolic-disease maps all embed the respiratory chain
and each lost exactly the same 11 oxidative-phosphorylation subunits (checked
directly, not inferred from the matching counts). The distinct signals are
glucosinolate biosynthesis (11.2x), bacterial chemotaxis (7.5x), beta-lactam
resistance (6.9x), ABC transporters (5.5x), 2-oxocarboxylic acid metabolism
(3.7x), oxidative phosphorylation (2.9x), two-component system (2.8x) and
secondary metabolite biosynthesis (1.7x) -- every one a large paralogous
family. Arbitration exists to fix over-calling caused by paralogy and
over-corrects hardest exactly where paralogy is highest: **it trades one
paralogy artefact for another.**

Per-gene gains, for the record, are individually significant: all seven
genomes have 95% CIs excluding zero for both precision and F1 (paired
bootstrap resampling genes, not calls), and precision improves 7/7
(one-sided sign test p = 0.0078).

On the three MAGs -- the actual target use case, where no ground truth exists
and none of this is verifiable -- arbitration resolves 2.07% of calls against
4.70% on reference genomes, i.e. its effect there is under half as large.
That bounds the untested exposure; it is not evidence it helps.

The mechanism is a level mismatch, not noise. A genome has many genes but one
KO set. A false-positive KO call on gene X usually names a KO genuinely
present on another gene, so dropping it does not shrink the set; dropping a
true positive that was a KO's only representative does. The precision gain is
roughly halved on the way up while the recall loss carries over intact.

Guidance added to `docs/methods.md`: use `best` when the per-gene assignment
is the product; **keep the default `all`** when the annotation feeds
presence/completeness analysis. No code change — the default was already
`all`, which is why this is a documentation correction and not a regression.
New `benchmarks/threshold_audit/downstream_impact.py` makes it reproducible.

## 3.21.0

**`--min-margin BITS`: a precision dial with a measured frontier.** Drops
calls closer than the given margin to their KO's own threshold. Combined with
`--multi-ko-policy best`, mean precision over the seven benchmark genomes:

| `--min-margin` | precision (no arb) | precision (**+ arb**) | recall | F1 |
|---|---|---|---|---|
| 0 | 0.874 | **0.903** | 0.886 | 0.893 |
| 10 | 0.894 | **0.915** | 0.864 | 0.888 |
| 20 | 0.907 | **0.924** | 0.840 | 0.879 |
| 50 | 0.932 | **0.943** | 0.739 | 0.826 |
| 80 | 0.947 | **0.956** | 0.630 | 0.755 |

Two results worth separating:

- **Arbitration is not the same thing as being stricter, and cannot be
  replaced by it.** Arbitration alone (precision 0.903, recall 0.886) beats a
  raised margin alone at 10 bits (0.894, 0.876) on *both* axes, and the two
  stack at every point on the frontier. Competing-KO error and weak-hit error
  are distinct failure modes.
- **F1 falls monotonically as the margin rises.** `--min-margin` is only worth
  using when a false positive genuinely costs more than a false negative,
  which is a property of the downstream analysis rather than of the tool. It
  is documented as a trade, not recommended as a better setting, and the
  default is 0 (off).

Verified end-to-end: `--multi-ko-policy best --min-margin 20` on
*M. jannaschii* gives precision 0.955 / recall 0.796, against KofamScan's
0.903 / 0.884 on the same input.

## 3.20.0

**`--multi-ko-policy best`: resolve competing KO calls, +2.9 points of
precision.** KOfam fits every KO's threshold **independently**, each
maximising its own F-measure; nothing in that procedure makes KOs compete. A
protein matching several related profiles clears all of their thresholds at
once, and KofamScan -- and therefore `mpph annotate`, which reproduces it --
emits every one.

KEGG's reference does not work that way. Across the seven benchmark genomes
spanning all three domains, **>=99.88% of genes carry exactly one KO** (max
observed 2), so the extras are over-calls almost by construction. They are
also where the errors are: on *E. coli* multi-KO genes are 13.5% of calls but
**69% of all false positives**, at precision 0.39 against 0.96 for single-KO
genes.

The new policy keeps, per gene, the KO furthest above *its own* threshold:

| genome | precision (all -> best) | F1 (all -> best) |
|---|---|---|
| *E. coli* K-12 | 0.900 -> **0.937** | 0.935 -> **0.941** |
| *B. subtilis* 168 | 0.858 -> **0.908** | 0.898 -> **0.911** |
| *M. jannaschii* | 0.903 -> **0.934** | 0.894 -> **0.904** |
| *Verrucomicrobia* sp. S94 | 0.845 -> **0.878** | 0.848 -> **0.861** |
| *Lentisphaerae* sp. WC36 | 0.848 -> **0.864** | 0.842 -> **0.849** |
| *S. cerevisiae* | 0.945 -> **0.960** | 0.911 -> **0.914** |
| *A. thaliana* | 0.815 -> **0.837** | 0.867 -> **0.874** |
| **mean** | **0.874 -> 0.903** | **0.885 -> 0.893** |

Precision *and* F1 improve on all seven, so this is a net gain rather than
recall traded away. On *E. coli* it puts `mpph annotate` ahead of the tool it
reimplements: F1 0.941 against KofamScan's 0.935.

- **Ranking is by margin above threshold, not raw score.** Thresholds span
  ~30 to >2000 bits, so raw scores are not comparable between KOs. Ranking by
  *relative* margin (score/threshold) was also tried and is consistently
  worse, contradicting the obvious scale argument -- so this was measured, not
  reasoned.
- **`--min-ko-gap BITS`** additionally requires the winner to beat the
  runner-up, dropping the gene when the evidence does not separate them. Mean
  precision reaches 0.909 at 40 bits with F1 flat, so it is a preference about
  error type rather than a better setting.
- **The default stays `all`.** Changing it would silently alter existing
  users' results and would invalidate the KofamScan-equivalence result. `best`
  is opt-in, and the two claims now stand together: identical to the incumbent
  by default, better than it on request.
- New `Assignment.margin` field records the distance above threshold, using
  whichever score the significance test actually compared (best-domain for
  `score_type = domain` KOs, full-sequence otherwise) -- `Assignment.score`
  remains the full-sequence score and for domain-scored KOs the two are
  genuinely different quantities.

Evidence and the comparison of five candidate rules:
`benchmarks/threshold_audit/`. The shipped implementation was validated
end-to-end against the offline analysis on three genomes, reproducing
precision, recall and F1 to six decimal places.

**Caveat that constrains the claim**: KEGG's one-KO-per-gene structure may be
partly a curation convention rather than pure biology, so some of this gain is
agreement with how the reference is built. What is not convention is KOfam
emitting several mutually exclusive orthology assignments for one protein --
that is a real consequence of fitting thresholds per KO in isolation.

## 3.19.1

**`mpph annotate`'s peak memory is set by `--cpus`, not by your input.**
3.19.0 established what the 2.3-11.8 GB peak is *not* (the target sequences)
but left what it *is* unexplained. Measured now, varying only the thread
count on one fixed 4,237-protein genome:

| `--cpus` | peak RSS | wall clock | speed-up |
|---|---|---|---|
| 1 | 1.0 GB | 21:07 | 1.0x |
| 4 | 2.9 GB | 5:17 | 4.0x |
| 28 | 11.8 GB | 2:07 | 10.0x |

Roughly 0.4 GB per thread, which accounts for the entire spread across the
benchmark without anything exotic. Two consequences worth acting on, now
documented in `--cpus`' own help text, `docs/`, and the benchmark README:

- **No proteome is too large on memory grounds**, and no large-memory node is
  needed. If RAM is tight, lower `--cpus`.
- **That costs little.** Parallel scaling is already well past its linear
  region at 28 threads (4 threads gives a perfect 4.0x speed-up; 28 gives
  only 10x), so 28 -> 4 threads cuts peak memory 4x for 2.5x the wall-clock.

**Correction to 3.19.0's own wording**: it reported prefetch 11.4 GB vs
stream 11.6 GB and called streaming "marginally worse". Two runs of the
*identical* prefetch command have since measured 11.4 and 11.8 GB, so that
0.2 GB gap is inside run-to-run noise and there is no difference to
attribute to sequence loading in either direction. The substantive claim it
supported -- that sequence loading is not what drives peak memory -- is
unchanged and now has a positive explanation rather than only a negative one.

Documentation and help text only; no change to any result.

## 3.19.0

**Validated beyond microbes: all three domains of life.** The toolkit was
written for, benchmarked on and documented around prokaryotes and MAGs.
Nothing in it is actually prokaryote-specific -- KO is a cross-domain
orthology system -- but "should work" is not a result, so the annotation
benchmark now includes *Saccharomyces cerevisiae* (6,021 proteins) and
*Arabidopsis thaliana* (48,265 proteins, 11x *E. coli*) against KEGG's own
KO assignments, under the same containers and the same KOfam database as
every other genome.

- **Accuracy does not degrade outside bacteria.** F1: yeast **0.911**,
  *Arabidopsis* **0.867**, against 0.894-0.934 for the prokaryotes.
  `mpph annotate` still reproduces KofamScan call for call -- Jaccard 1.000
  on yeast and 0.9999 on *Arabidopsis*.
- **The residual disagreement is fully explained.** All three differing calls
  across *Arabidopsis*' 12,223 lie within 0.016 bits of the KO's threshold,
  in both directions. HMMER prints bit scores to one decimal and KofamScan
  parses that text, so it compares a rounded score against a two-decimal
  threshold; `mpph annotate` compares pyhmmer's in-memory value at full
  precision. The two agree on every call that is not a tie inside the
  reference database's own printing precision.
- **Coverage does fall, and that is biology**: 76% of *E. coli* proteins get
  a KO, 62% of yeast, 24% of *Arabidopsis* -- for all three tools alike. KO
  describes metabolism and core cellular processes, not a plant's whole
  protein complement. Recall against what KEGG *does* assign is the highest
  in the benchmark (0.925 for *Arabidopsis*).
- **Eukaryotic isoforms are now scored fairly.** A eukaryotic proteome lists
  every splice variant while KEGG names one representative protein per gene.
  Without a restriction to the evaluable set, 10,025 calls on isoforms the
  reference cannot represent would have counted as false positives, dragging
  a hypothetical perfect tool to precision 0.5. Applied identically to all
  three tools; a no-op for prokaryotes.
- README, `docs/`, and `paper/manuscript.md` rewritten from "microbial /
  MAG" framing to any genome, with the eukaryote-specific caveats stated
  where a user would hit them. Figure 1 gains the two eukaryotes and a
  domain-of-life column.

**Correction to 3.18.x's `--sequence-loading` rationale.** That flag shipped
with the claim that pyhmmer's prefetching drives peak memory and would not
scale to a eukaryotic proteome. Measured, that is wrong in both halves:

- A prefetched sequence block costs **~1.0 kB per protein** -- 4.5 MB for
  *B. subtilis*, 45.5 MB for all of *Arabidopsis*. The earlier reasoning
  extrapolated a roughly constant quantity as if it scaled with input size.
- Re-running *B. subtilis* both ways: prefetch 11.4 GB / 2:06.9, stream
  11.6 GB / 2:07.5, byte-identical output. Streaming is marginally *worse*.

The flag is kept -- it is a real control for metagenome-scale protein
catalogues, where ~1 kB each does add up -- but `auto` now switches at
**1,000,000** proteins rather than 10,000, and the documentation says what
was measured instead of what was assumed. Peak memory (2.3-11.8 GB) is
dominated by the profile search, not the targets: the worst case in the
whole benchmark is the 4,237-protein *B. subtilis* at 11.8 GB, while the
48,265-protein *Arabidopsis* peaks lower at 9.7 GB.

## 3.18.1

Documentation and validation-figure release; no package code changed.

- New `paper/manuscript.md`: an Applications Note draft describing what this
  tool is now. The 2023 preprint's framing ("constructing phylogenetic trees
  based on metabolic pathways") no longer matches a 13-subcommand toolkit
  that annotates, profiles, tests and draws. Every figure in it is
  reproducible from `benchmarks/`; the outstanding decisions (authorship,
  target venue, whether to add a PVC application section) are listed
  explicitly in the draft rather than silently assumed.
- New `benchmarks/figures/fig3_phylo_calibration`: Type I error uncorrected
  vs corrected, one point per simulated no-effect scenario.
- **Correction to 3.18.0's reported numbers**: `benchmarks/phylo/README.md`
  said 41 of 50 cells exceeded the nominal 0.05 under the uncorrected test.
  The real count is **31 of 50** (verified against `calibration.csv`). Every
  other figure in that table -- worst 0.948 vs 0.072, median 0.086 vs 0.000,
  1 corrected cell above nominal -- was correct as published.

## 3.18.0

**`mpph compare --tree`: phylogeny-aware testing.** `compare`'s Fisher and
Mann-Whitney tests assumed every genome is an independent observation. They
are not -- close relatives share features by common descent -- so the tests
were anti-conservative, which this package had documented as a caveat
("treat as exploratory") in four places rather than fixed. Measured across a
54-cell grid of tree shapes, prevalences and evolutionary rates with **no
group effect present at all**, the uncorrected test's Type I error reaches
**94.8%** at a nominal 5%.

- New `mpph.phylo`: a Newick parser that keeps branch lengths (the existing
  `treecompare` one is topology-only by design), tip pruning that folds a
  collapsed degree-2 node's length into its surviving child, Brownian and
  symmetric-Mk simulation, Felsenstein pruning with per-node rescaling, ML
  rate fitting, and Fitch parsimony.
- `differential_features(..., tree=...)` and `mpph compare --tree FILE
  --phylo-permutations N --seed N` add `p_value_phylo`/`q_value_phylo`,
  `n_state_changes`, `phylo_confounded`, `n_phylo_replicates_accepted`,
  `phylo_k_tolerance` and `phylo_warning`. **Without `--tree` nothing
  changes**, columns included.
- **The presence null is conditioned on the observed number of present
  genomes.** This is what makes the test work rather than a refinement:
  unconditioned, about a third of simulated replicates come out invariant,
  contribute a zero statistic, dilute the tail, and the test then calls a
  purely clade-confounded feature significant (p = 0.013 at 16 tips, getting
  worse as the tree grows -- i.e. failing on exactly the case it exists to
  catch). Conditioned, the same input gives p ~ 0.78.
- **The completeness null fits nothing.** Scored with a rank statistic
  (Cliff's delta), which is invariant to positive scaling and translation,
  the Brownian null depends on neither the rate nor the ancestral state --
  so no optimizer, and one null serves every feature.
- **mpph's own dendrogram is rejected as `--tree`**, by topology rather than
  filename so a renamed copy is caught too: it is built from the very
  features being tested, so correcting those tests with it is circular.
  A cladogram without branch lengths is also rejected rather than defaulted
  to 1.0, and a tree matching under 80% of the organisms errors instead of
  silently analysing the survivors.
- `--phylo-permutations` defaults to **9999, not 999**: a permutation q
  cannot fall below `n_features/(replicates+1)`, so at 999 a 400-feature
  matrix could never produce a q under 0.40 -- a purely arithmetic artifact
  that reads as "nothing survives correction". `compare` warns when the
  floor still exceeds 0.05.
- `n_state_changes` (Fitch, exact and simulation-free) is the number that
  explains the result: when it is 1 the feature arose once, and if that
  origin sits on the branch separating the groups, no method can separate
  association from coincidence (Maddison & FitzJohn 2015). A large
  `p_value_phylo` there is the correct answer.
- Calibration gate in `benchmarks/phylo/`: worst corrected Type I 0.072
  against 0.948 uncorrected, median 0.000 against 0.086, with 1 of 50 cells
  above nominal where 3.9 are expected by chance. Regression guards live in
  `tests/test_phylo.py`, including hand-computed pruning likelihoods
  (0.2161661792 / 0.2838338208), a bit-identical-across-rootings check, an
  exact VCV, and the clade-confounded acceptance test.

## 3.17.1

New `benchmarks/` directory: head-to-head validation against the established
tools this package overlaps with, on real data. **No package code changed in
this release** -- the benchmarks were run against 3.17.0 and apply unchanged
to 3.17.1; only `benchmarks/` and a `.gitignore` rule are new.

- **Annotation** (`benchmarks/annotation/`) -- `mpph annotate` vs KofamScan
  1.3.0 vs eggNOG-mapper 2.1.15, over three difficulty tiers: three model
  organisms (*E. coli*, *B. subtilis*, *M. jannaschii*), two non-model
  environmental genomes KEGG never hand-curated (*Verrucomicrobia* sp. S94,
  *Lentisphaerae* sp. WC36), and three real marine MAGs with no reference
  annotation at all. `mpph annotate` reproduces KofamScan at Jaccard 1.000
  on six of the eight genomes (0.9997 and 0.9994 on the other two) while
  running against the same KOfam database; F1 against KEGG's own KO
  assignments is 0.842-0.934. The profile-HMM tools' margin over
  eggNOG-mapper is *larger* on the non-model genomes than on the model
  organisms.
- **Enrichment** (`benchmarks/enrichment/`) -- `mpph gsea`/`mpph enrich` vs
  clusterProfiler 4.18.4 on identical `TERM2GENE` category definitions and
  identical size filters: GSEA enrichment-score Spearman rho = 1.000
  (NES 0.999, p-value 0.994), ORA p-value rho = 1.000.
- **Reported against this package, not omitted**: `mpph annotate`'s peak
  memory is 10-17x KofamScan's (2.3-12.3 GB vs 0.13-0.72 GB) and scales
  erratically with input size. See `benchmarks/annotation/README.md`.
- Both comparison harnesses ship with tests that check their metrics against
  hand-constructed inputs with known expected values (`pytest benchmarks/`,
  deliberately outside the package's own `testpaths`). Figures regenerate
  from the committed summary tables via `python benchmarks/make_figures.py`,
  with no reference database needed.
- `.gitignore`: the repo-wide `*.csv` rule was hiding the benchmark summary
  tables, which are the evidence behind the README claims and the input to
  the figure script -- exempted.

## 3.17.0

**`--top-category`** for `mpph enrich`/`mpph gsea` (`--ontology kegg-pathway`
only): restrict the tested universe to one BRITE top-level category (e.g.
`Metabolism`). KEGG PATHWAY isn't metabolism-only -- it also covers Genetic
Information Processing, Environmental Information Processing, Cellular
Processes, Organismal Systems, Human Diseases and Drug Development, so an
unrestricted run can surface a significant hit from any of these. `mpph run`
already had this exact filter (`--top-category`/`--all-categories`) for the
presence heatmap; `enrich`/`gsea` didn't.

- Reuses the same `mpph.pathways.fetch_pathway_categories()` (`br08901`)
  source `run --top-category` already uses -- no new KEGG endpoint.
- Unset by default: **no behaviour change** for existing usage -- every
  category is still tested unless `--top-category` is explicitly given
  (deliberately not mirroring `run`'s own default-to-Metabolism behaviour,
  since `enrich`/`gsea` are already-shipped subcommands and changing their
  default would silently change existing results).
- The output CSV always gains a `top_category` column when `--ontology
  kegg-pathway`, whether or not the flag is used, so the current, unrestricted
  default output can still be filtered by category afterward without rerunning.
- A `--top-category` value with `--ontology kegg-module`/`go` (which have no
  BRITE top-level category) prints a warning and is otherwise ignored, rather
  than silently doing nothing.

## 3.16.3

Lint cleanup + ruff upgrade: fixed the 44 pre-existing style findings that
3.16.1 pinned around instead of fixing, then removed the `<0.16` upper
bound so `dev`'s `ruff` dependency tracks the linter's latest release
again. Each finding was checked individually rather than run through a
blanket `--fix`, because ruff's isort (`I001`) and dunder-all (`RUF022`)
autofixes both collapse comment-delimited groupings into one
alphabetically-sorted block. `mpph/__init__.py` uses exactly that pattern
for its KEGG/data-access/enrichment/gsea/kgml/matrix+plot import sections
and its `__all__` category comments -- a blind `--fix` would have
scattered unrelated symbols under the wrong header (e.g. `KGMLNode`
sorting in under the "KEGG connection interface" comment meant for
`kegg.py`'s exports). Import sections now use `# isort: split` between
groups, which ruff honors for independent per-section sorting without
touching the sections' relative order; `__all__` keeps its category
comments behind a scoped `# noqa: RUF022`, since ruff's own docs describe
preserving comment-delimited categories as out of scope for that rule's
fix. Also fixed: a genuine duplicate `.modules` import in `cli.py`
(`I001`); several `# noqa` suppressions for bandit/pycodestyle codes
(`S314`, `S202`, `S307`, `E402`) that this project's ruff config doesn't
enable and were therefore dead (`RUF100`); two `sub == sub` /
`score != score` NaN self-comparison idioms replaced with `math.isnan()`
(`PLR0124`); a `ValueError` -> `TypeError` fix in `traits.py`'s
panel-shape check (`TRY004`), which also required adding `TypeError` to
`cli.py`'s top-level error handler so a malformed panel file still prints
a clean `Error: ...` message instead of an unhandled traceback; and, in
the test suite, several mock `Session` classes' `headers: dict = {}`
class attributes (`RUF012`) moved into `__init__` rather than annotated
`ClassVar`, since `mpph.kegg.kegg_get` really does call
`session.headers.update(...)` and a shared class-level dict would have
leaked mutated headers across test instances.

## 3.16.2

CI hotfix #2: three of `tests/test_cli_annotate.py`'s tests
(`..._missing_kofam_db_gives_clean_error_not_traceback`,
`..._happy_path_writes_output_and_manifest`,
`..._default_out_path_is_fasta_stem_annotated`) passed locally (where
`pyhmmer` happens to be installed) but failed in CI, where the `annotate`
extra is deliberately not part of the `dev` install: `cmd_annotate`'s own
`import pyhmmer` availability probe runs before any of the code paths these
tests exercise, so a genuinely-missing `pyhmmer` short-circuits to the
"install the extra" message before ever reaching them. `tests/
test_kofam_pyhmmer.py` already had the right guard for this
(`pytest.importorskip("pyhmmer")`); these three needed the same guard and
didn't get it. Fixed by adding it to each; verified by simulating a
missing `pyhmmer` locally (blocking the import) and confirming these three
now skip instead of fail, while the four tests that don't need real
`pyhmmer` (argparse wiring, `--setup-db`, the missing-dependency message
itself) still run and pass either way.

## 3.16.1

CI hotfix, unrelated to 3.16.0's actual content: `dev`'s `ruff>=0.1` had no
upper bound, so CI's fresh `pip install -e ".[dev]"` picked up ruff 0.16.0
the moment it was released and started failing on 44 pre-existing style
findings across the repo (comment-divided import groups, an unused unpacked
variable, a `dict.items()` that only uses the values, etc.) that every
prior ruff 0.15.x release never flagged -- none introduced by 3.16.0's own
changes (verified: ruff 0.15.22, the latest patch below the new release,
passes the repo cleanly both before and after 3.16.0's diff). Pinned
`ruff>=0.1,<0.16` so a future linter release can't silently break CI again
without a deliberate version bump. The 44 pre-existing findings themselves
are unaddressed -- a separate cleanup, not bundled into this hotfix or into
3.16.0's `mpph annotate` work.

## 3.16.0

New subcommand: **`mpph annotate`** -- locally annotates a protein FASTA
with KEGG KO ids, entirely offline after a one-time database download. No
external HMMER/KofamScan/eggNOG-mapper install needed.

`mpph` previously assumed annotation already happened elsewhere (`--codes`
for a KEGG organism, `--user` for your own already-annotated KOs). This
closes that gap the way KOfam itself generalizes -- across the tree of
life, not a fixed set of reference species (the approach tools like KOBAS
use, which degrades with evolutionary distance from whichever reference is
closest).

- New `mpph.kofam` module: `annotate_fasta` searches every KOfam HMM
  profile against every input protein via
  [pyhmmer](https://pyhmmer.readthedocs.io/) (Cython bindings to HMMER3,
  `pip install pyhmmer` -- no separate HMMER binary/compilation needed),
  reimplementing KofamScan's own scoring/assignment logic rather than
  shelling out to compiled HMMER + KofamScan's Ruby wrapper.
- Significance rule ported from KofamScan's actual source
  (`result/hit.rb`), not its README, which disagrees with its own code on
  one point: the comparison is `score >= threshold`, not "higher than".
  `score_type == "domain"` KOs compare their best single domain's score;
  every other KO compares the full-sequence score (`result/parser.rb`). KOs
  with no threshold in `ko_list` (too few reference sequences) are never
  assignable, by KOfam's own design.
- `download_kofam_db()` / `mpph annotate --setup-db DIR`: downloads
  `ko_list.gz` + `profiles.tar.gz` from KEGG's own KOfam distribution
  (~1.5 GB compressed) over HTTPS, skips re-downloading if `DIR` already
  looks populated.
- Output is `gene<TAB>K#####`, one row per significant assignment, no
  header -- the same shape as KofamScan `--format mapper`, so it's already
  readable by `--user` with zero adapter code (verified end-to-end against
  the real, unmodified `load_user_kos()`, not just asserted by shape).
- New `annotate` extras group (`pip install mpph[annotate]`) -- `pyhmmer`
  is deliberately not part of the base install or the `dev` extra, so
  `pip install mpph` stays thin; running `annotate` without it prints an
  actionable error instead of a traceback.
- Scope for this release is **KO-only, no GO** (GO would need eggNOG's
  differently-structured per-clade databases, a separate addition) and
  **protein FASTA only, no ORF prediction** (run Prodigal or similar first
  for raw contigs).

## 3.15.0

**`--unknown-policy`**: `mpph compare` and `mpph pan` now let you choose how
`NaN` ("unknown", never "confirmed absent") is treated by the actual
statistics, instead of always silently excluding it. This is distinct from
the 3.10.1 heatmap hotfix, which only fixed how `NaN` is *rendered* -- the
computations behind `differential_features`/`pan_classify` still hardcoded
"drop it" with no way to opt into "count it as absent" or demand "refuse to
run" instead. Completes the last unaddressed item from the third external
review.

- `differential_features()`/`pan_classify()` gain `unknown_policy`
  (`"exclude"` default / `"absent"` / `"error"`):
  - `"exclude"`: known-only denominator -- unchanged default behaviour.
  - `"absent"`: explicit opt-in to counting `NaN` as absent. Never the
    default, since it may misrepresent an assembly/annotation gap (e.g. an
    incomplete MAG) as genuine absence -- only use it with a specific
    reason to believe that's actually true here.
  - `"error"`: raise immediately if any `NaN` is present among the compared
    groups/matrix, forcing missingness to be resolved (e.g. via
    `--min-genome-completeness` QC filtering) before you compare or
    classify at all.
  - `n_known`/`n_unknown`/`known_fraction` (and the min-known thresholds
    that gate `insufficient_known_values`/`insufficient-data`) are always
    computed from the **true, policy-independent** evidence -- a feature
    with too little real data is never waved through just because
    `"absent"` would otherwise happily fill the gaps with zeros.
- `pan_classify()` gains a new `n_absent` output column (the complement of
  `n_present` under whichever denominator `unknown_policy` selects) and a
  `min_known_samples` threshold alongside the existing `min_known_fraction`
  (both independent of `unknown_policy`, for the same reason as above).
- New CLI flags: `mpph compare --unknown-policy {exclude,absent,error}
  --min-known-samples N`; `mpph pan --unknown-policy {exclude,absent,error}
  --min-known-fraction F --min-known-samples N` (`--min-known-fraction` was
  already a `pan_classify()` parameter but had no CLI flag until now).
- Incidental fix: regenerating `examples/Synechococcus_pan_classes.csv` for
  the new `n_absent` column surfaced a pre-existing inconsistency -- its
  `feature_id`s (`M00019`, `M00049`, ...) were KEGG module identifiers that
  could only have come from a completeness/module-mode matrix, while the
  currently-committed `Synechococcus_matrix.csv` (and its manifest's
  `"mode": "presence"`) is pathway-based. No module-mode source data exists
  anywhere in the repository to reproduce the original file's exact
  content, so it is regenerated directly from the current
  `Synechococcus_matrix.csv` instead -- making the file internally
  consistent with its own committed source for the first time.

## 3.14.0

New subcommand: **`mpph pathmap`** — draws a KEGG pathway (e.g. `ko00010`,
Glycolysis) or the entire global metabolic network map (`ko01100`, every
KEGG pathway stitched into one diagram, ~3,800 reactions) using KEGG's own
KGML layout, colouring each enzyme node and the reaction(s) it catalyzes by
whether it's present in group A, group B, both, or neither — comparative,
not KEGG's own static colouring.

- New `mpph.kgml` module: `fetch_kgml`/`parse_kgml` (KGML XML -> positioned
  compound/enzyme/map-reference nodes + the reactions linking them),
  `normalize_map_id` (`"01100"`/`"map01100"`/`"path:ko01100"` all resolve to
  the `ko#####` form KEGG actually serves KGML for), `reactions_for_ortholog`.
- New `plot_kgml_map`: an enzyme entry counts as present in a group if *any*
  of its KOs are (KEGG sometimes lists several isozymes under one entry);
  a reaction catalyzed by more than one such entry takes the union of all
  their KOs. Every substrate-product pair of a reaction is drawn (not just a
  guessed "primary chain"), since KGML's own compound order isn't a reliable
  signal and some real steps genuinely branch (e.g. aldolase, one substrate
  to two products). Compound nodes and map-reference boxes are drawn from
  KGML's fixed layout for context, not coloured by data.
- `--codes-a`/`--codes-b` (organism codes, KOs unioned) or `--user-a`/
  `--user-b` (your own annotations) per group; a manifest records both
  groups' sources/KO counts, per-status enzyme counts, and the usual
  provenance fields.
- New example: `examples/Prochlorococcus_MIT9313_vs_AS9601_pathmap.png` —
  the same low-light/high-light strain pair as the existing enrichment
  example, now as a network diagram. Porphyrin (chlorophyll/heme)
  metabolism comes back entirely shared between the two strains — the core
  photosynthetic-pigment pathway both need, distinct from the antenna-protein
  genes the enrichment example found unique to the low-light strain.

## 3.13.0

Manifest reproducibility: `run` and `traits` now record the calling
directory's git commit SHA (`git_commit`, `null` outside a checkout -- a
best-effort lookup that never fails the run) and the core dependency
versions actually used (`dependency_versions`: numpy, pandas, scipy,
requests). Together with the already-recorded `mpph_version`/`python`/
`kegg_release`/full `command`, a manifest now identifies exactly which
version of everything produced a given result.

This completes the engineering batch (cache collision-proofing + atomic
writes, HTML report large-matrix protection, `cliffs_delta` performance,
manifest provenance) from the third external review's remaining
recommendations, on top of the P0 fixes (3.5.1-3.5.3) and P1 scientific-
rigor batch (3.6.0-3.10.1). Code reorganization into subpackages was
explicitly not part of this batch. The release/repo-management items
(GitHub Release, wheel-based CI, Docker, Dependabot) remain untouched and
need explicit sign-off before any of that is attempted.

## 3.12.1

`cliffs_delta` performance: replaced the O(n*m) pairwise comparison matrix
(`a[:, None] > b[None, :]`) with the algebraically equivalent Mann-Whitney U
statistic (`delta = 2U/(n*m) - 1`, computed via sorting in O(n log n + m log
m)) -- for two groups of 20,000 each, ~95x faster and ~400MB less peak
memory in this synthetic benchmark. Matters once a group is thousands of
features/genes rather than the usual handful of organisms. Verified
mathematically identical to the direct pairwise count, ties included (not
just empirically checked -- `2U/(n*m) - 1` and `(greater - less)/(n*m)` are
the same quantity, since ties count as 0.5 in U the same way they count as
neither "greater" nor "less" directly). The one shipped example with a
`cliffs_delta` column did not need regenerating, since the function's
output for the same input is unchanged.

## 3.12.0

HTML report: large-matrix protection. The interactive grid embedded the
full `organisms x features` matrix as JSON and built one styled DOM node per
cell with no size limit at all -- for a large analysis (thousands of
organisms and/or features), that's hundreds of thousands to millions of DOM
nodes, easily enough to hang or crash the browser tab.

- New `build_report(..., max_cells=50_000)` / `mpph report --max-cells N`:
  above this many cells, the matrix is not embedded in the JSON payload at
  all (keeping the HTML file itself small) and the grid section shows a
  plain-text summary (dimensions, cell count, pointers to the CSV/figure
  outputs) instead of attempting the table.
- The QC and manifest tabs are unaffected either way -- they're small
  regardless of matrix size -- and the full data is always in
  `<slug>_matrix.csv` / `_ordered_matrix.csv` and the heatmap figure(s).
- Every existing shipped example is well under the default threshold, so
  nothing needed regenerating; this closes a real gap for any larger
  analysis, including a user's own `build_report()` call.

## 3.11.0

KEGG response cache: collision-proof filenames and atomic writes. Reproduced
with a failing case first, now has regression tests. **Upgrading invalidates
any existing `.mpph_cache/` directory** (filenames changed) -- harmless,
since the cache is a pure performance optimization; the next run just
re-fetches from KEGG.

- `_cache_path` collapsed every run of non-alphanumeric characters in an
  endpoint to a single `_`, so `"a/b"`, `"a//b"` and `"a_b"` (a literal
  underscore is non-alphanumeric too) all produced the identical cache
  filename -- a future endpoint construction that happened to collapse the
  same way would silently serve one query's cached response for a different
  one. Now a readable slug plus a hash of the untransformed endpoint string,
  which can't collide regardless of what the slug collapses to.
- Cache writes went straight to the final filename; a process killed
  mid-write (or two processes racing the same cache entry) could leave a
  truncated file that a later read would silently treat as a complete KEGG
  response. Now written to a per-process temp file first, then moved into
  place with `os.replace()` (atomic on both POSIX and Windows).
- The 11 committed offline-test cache fixtures are renamed to match (same
  content, new filenames); the offline test suite (`test_integration.py`,
  `test_cli_enrich_offline.py`, `test_cli_gsea_offline.py`) still passes
  entirely from cache with no live network call.

## 3.10.1

Heatmap rendering hotfix: `NaN` ("unknown / not assessed") is now visually
distinct from both "confirmed absent" and every real value, in both modes.
Reproduced with a failing case first, now has regression tests.

- `mode='presence'`: a `NaN` cell was indistinguishable from a confirmed-
  absent one (`NaN > 0` evaluates to `False`, the same as `0 > 0`), so an
  unknown pathway silently looked exactly like a confirmed-absent one.
- `mode='completeness'`: an undetermined module score (the 3.6.0 tri-state
  scoring fix means `module_completeness` can now genuinely return `NaN`)
  rendered as **solid black** -- the colormap correctly marks `NaN` with a
  fully transparent "bad" colour, but the code discarded the alpha channel
  before drawing, leaving opaque black behind. This was worse than looking
  "absent": a stark, attention-grabbing artifact with no legend explanation.
- Both now render in a new, distinct colour (`UNKNOWN`, a muted lavender-gray
  that doesn't collide with any category hue or with `ABSENT`), with a
  legend entry ("unknown / not assessed") that appears only when the plotted
  matrix actually contains a `NaN`.
- No shipped example needed regenerating: none of the committed example
  matrices contain `NaN` today (`build_presence_matrix` never produces it,
  and none of the completeness-mode examples currently hit an unresolved
  module reference) -- this closes a real gap for any matrix that does,
  including a user's own data passed to the public `plot_matrix` API.

## 3.10.0

MAG genome-quality metadata wired into the `run` pipeline. `mpph.samplesheet`
already had `import_checkm2()`/`import_gtdbtk()` parsers, but nothing in the
CLI ever called them -- a missing feature in an incomplete MAG could look
like true absence with no way to filter, annotate, or even be warned about it
from the command line.

- New `run --qc-metadata FILE`: a table keyed by `sample_id`/organism with a
  `completeness` column (optionally `contamination`, `taxonomy` -- e.g.
  `import_checkm2()`/`import_gtdbtk()` output saved to TSV). Joined into
  `*_qc.csv` and recorded in the manifest.
- New `--min-genome-completeness N`: drops organisms below `N`% completeness
  before scoring, reported the same way as `--keep-empty`'s exclusions.
  Any *included* organism under 90% completeness is warned about regardless.
- New `mpph.samplesheet.mimag_quality_tier()`: the standard MIMAG (Bowers et
  al. 2017) high/medium/low tiers from completeness+contamination (the full
  standard also needs rRNA/tRNA evidence this doesn't assess -- documented
  as such, not overclaimed). `import_checkm2()` now includes a
  `quality_tier` column using it.
- Documented explicitly, in both `docs/methods.md` and the new
  `apply_genome_completeness_qc()`'s docstring: this filters or flags, it
  does **not** attempt a linear completeness correction
  (`observed / genome_completeness`) -- different functional genes are not
  lost uniformly at random as assembly quality drops, so that correction
  would fabricate precision the data doesn't support.

This completes the P1 scientific-rigor batch (Module AST, GSEA, ORA, trait
panels, tree/RF, MAG QC) from the third external review's most severe
correctness/rigor items, on top of the three P0 fixes in 3.5.1-3.5.3. The
review's remaining ~11 items (absent/unknown heatmap semantics, manifest/
cache/report engineering, code reorg into subpackages, `cliffs_delta`
performance, and release/repo-management asks) are still unaddressed and
were not part of what was asked for in this round.

## 3.9.0

Newick label safety and richer tree comparison. Reproduced with a failing
case first, now has regression tests.

### Newick export could silently collide two different labels
- `linkage_to_newick`'s label sanitizer deleted/replaced characters
  (`(),:;[]'` stripped, space -> `_`) instead of using Newick's own quoting
  rule -- so two genuinely different labels could produce the identical
  output name. A real, not just hypothetical, case: the microbiology
  convention `"[Eubacterium] rectale"` (bracketed provisional genus) and a
  plain `"Eubacterium rectale"` both became `"Eubacterium_rectale"`. Also,
  every organism label this tool itself produces (`"Name (code)"`) contains
  parens, so this affected the tool's own default output, not just edge
  cases.
- Fixed to follow the Newick spec: a label needing special handling is
  wrapped in single quotes (internal `'` doubled) and left verbatim,
  preserving it losslessly; a plain space with nothing else special still
  uses the conventional (unquoted) underscore substitution for the common
  case. `linkage_to_newick` now also raises `ValueError` naming the two
  labels if a collision survives even this (the one case quoting alone
  can't disambiguate: a literal underscore vs. a space-turned-underscore) --
  refusing to silently merge two samples' identities in the exported tree.
- `treecompare`'s Newick parser is now quote-aware (handles a quoted token,
  including the doubled-quote escape), so it correctly reads trees this
  package itself exports as well as reference trees from other tools that
  use standard Newick quoting.
- Both example organism trees (`Prochlorococcus_organism_tree.nwk`,
  `Synechococcus_organism_tree.nwk`) are regenerated with the new quoting
  (same topology and branch lengths, only the label formatting changed) --
  one of the regenerated labels, `Synechococcus sp. JA-2-3B'a(2-13) (cyb)`,
  already contains a literal apostrophe, which the old code silently
  dropped rather than escaping.

### Richer, more explicit `robinson_foulds` output
- Added `only_in_a`/`only_in_b` (the actual leaf-name lists, alongside the
  existing `n_only_in_a`/`n_only_in_b` counts) and `rooted: true`, making
  explicit that clades are compared as rooted descendant-sets (appropriate
  for a UPGMA dendrogram and a typically-rooted reference tree) rather than
  unrooted bipartitions-up-to-complement.

## 3.8.0

Trait panel versioning and provenance. `mpph traits` had no manifest at all
(unlike `run`, `enrich` and `gsea`), so a trait-scoring run left no record of
which panel file, version, or KEGG release produced it.

- Built-in panels (`biogeochemistry`, `respiration`, `carbon_fixation`) now
  declare `panel_version` (starting at `1.0.0`) alongside the existing
  `schema_version` and the "illustrative, verify before publication"
  disclaimer.
- New `panel_provenance()` returns the resolved panel path, its SHA-256
  content hash, and whatever top-level metadata fields the panel file
  declares -- a hash lets a later run detect that a panel's content changed
  underneath a previously-recorded analysis. `load_trait_panel`'s existing
  signature and behavior are unchanged (refactored to share path/file
  resolution, not rewritten).
- `mpph traits` now writes `<slug>_manifest.json` (command, KEGG release or
  `n/a (user data)`, organism/trait counts, panel provenance, outputs),
  matching the manifest already written by `run`/`enrich`/`gsea`.

This does **not** implement the fuller per-marker provenance schema some
reviews suggest (curator names, literature DOIs per marker, taxonomic scope,
paralog caveats) -- that content doesn't exist anywhere in this repository
today, and fabricating placeholder citations would be worse than not having
them. Filling that in for the 3 built-in panels is future curatorial work,
not something to synthesize.

## 3.7.0

ORA (`mpph enrich`) correctness fix and richer output. Reproduced with a
failing case first, now has regression tests.

### `ko:` namespace prefix silently caused zero matches
- A study/background list copied straight from KEGG's own `link/ko/<org>`
  endpoint, or an eggNOG-mapper `KEGG_ko` column, keeps the `ko:` prefix
  (`ko:K00001`) -- but this module's own category-membership fetchers key on
  the bare id (`K00001`). The mismatch wasn't an error, just a silent
  all-unannotated result (0 categories tested) for anyone whose input
  happened to use that natural, KEGG-native format. Both `study` and
  `background` are now normalized (`ko:K#####` -> `K#####`) before matching;
  GO-mode gene ids are untouched (only that exact pattern is stripped).

### New: `odds_ratio` / `odds_ratio_ci_low` / `odds_ratio_ci_high`
- The standard 2x2-table effect size (Haldane-Anscombe corrected for a zero
  cell), more comparable across categories of very different size than
  `fold_enrichment` alone (a ratio of rates, which can look arbitrarily large
  for a tiny category with a single lucky hit). Added to
  `examples/Prochlorococcus_MIT9313_unique_enrichment.csv`, computed directly
  from that file's own already-published counts (no re-run needed, so every
  other value is unchanged).

### Documentation
- Both the module docstring and `docs/methods.md` now state plainly that GO
  enrichment here is **flat** (no GO-DAG ancestor propagation) -- a gene
  counts only toward the terms your mapping lists for it, not their
  ancestors, unless the mapping already includes them. This is a difference
  in method from tools like topGO/clusterProfiler-with-`.db`, not a bug.

## 3.6.1

GSEA correctness and memory. Reproduced with a failing case first, now has
regression tests.

### Non-finite scores rejected
- `load_ranked_list`: `float()` accepts `"nan"`/`"inf"`/`"-inf"` without
  raising, so these previously entered the ranking silently -- an `Inf`
  score sorts to one end and dominates the running-sum statistic. A
  realistic real-world source: DESeq2 and similar tools report `Inf`/`-Inf`
  log-fold-change or `NA` for genes with zero counts in one group, and this
  tool's own docs recommend piping such a tool's statistic in directly via
  `--ranked-list`. Now rejected with a clear error naming the offending
  genes. The bundled `examples/Ecoli_ranked_logFC.tsv` itself had 25 such
  genes (of 3047) from the source study's own logFC column; regenerated with
  them dropped (`examples/Ecoli_ranked_logFC_gsea.csv` updated to match --
  the top hit, Ribosome, is unchanged at NES=-1.96, q=0.0054).
- `rank_from_expression`: `dropna()` doesn't catch `+-Inf` surviving the
  epsilon-padding meant to prevent literal division by zero (e.g. an `Inf`
  already present in the input expression matrix). Now rejected the same way.

### Memory
- The gene-set permutation null distribution allocated one
  `(permutations, n_genes)` array (plus several same-shaped temporaries) per
  distinct category size -- for tens of thousands of genes and thousands of
  permutations this can reach multiple GB at once. Now processed in batches
  of 200 permutations, bounding peak memory regardless of the total
  permutation count. Verified to produce bit-identical results to the
  unbatched form for the same seed (numpy's `Generator` is stream-based, so
  requesting the same draws in smaller chunks doesn't change them).

### Documentation
- Sharpened the module docstring: this is gene-set permutation on a fixed
  ("preranked") ranking giving a *nominal* p-value, not the original GSEA's
  default *phenotype* permutation (which re-derives the ranking from the raw
  samples and so also captures gene-gene correlation structure) -- phenotype
  permutation is not implemented here. `docs/methods.md` already documented
  this distinction from an earlier pass.

## 3.6.0

Module-completeness scoring now distinguishes "confirmed absent" from
"can't tell" -- a behavior change (affected modules now score `NaN` where
they previously scored a plausible-looking but fabricated number), so this
is a minor bump rather than another patch. Reproduced with a failing case
first, now has regression tests.

### Module completeness: unresolved references are unknown, not absent
- A step whose truth value depends on a nested `M#####` reference with no
  definition available -- not fetched, filtered out by the `Pathway`-only
  default (`--all-modules` off), or a cyclic reference -- silently counted
  as "not satisfied" (0), pulling completeness down as if the step were
  confirmed missing. Same for a step this parser's expression evaluator
  can't handle. Both are now correctly treated as **unknown**, excluded from
  both the numerator and denominator (the same treatment `--` placeholder
  steps already got), using two-valued (Kleene) evaluation: an undetermined
  token is substituted with both 0 and 1, and only a genuine disagreement
  (the step's truth value actually depends on it, e.g. an AND with no other
  determined factor) is reported as unknown -- an OR already satisfied by a
  known KO, or an AND already falsified by a known-missing KO, is still
  correctly determined regardless of the unresolved reference.
- `module_completeness` now returns `NaN` (never `0.0`) when every real step
  in a module turns out undetermined -- previously a purely circular or
  entirely-unresolvable definition silently scored 0.0, indistinguishable
  from genuine absence.
- `evaluate_module`'s `state` is now `"unknown"` (not derived from a score
  that already baked in "unresolved = failed") when nothing is determinable;
  `StepResult.satisfied` is `bool | None`, and `mpph explain` marks an
  unknown step `??` instead of showing it as `--` (confirmed missing).
- `filter_matrix` (the `--min-prevalence`/`--max-prevalence`/`--drop-core`
  machinery) now computes prevalence among *known* organisms only, so a
  module that is undetermined for some organisms isn't silently counted as
  absent in them; a module undetermined for *every* organism has no
  computable prevalence and is dropped rather than crashing.

## 3.5.3

Correctness hotfix for user-input loading. Reproduced with a failing case
first, now has a regression test.

### `--user` input auto-detection
- A single-genome KofamScan-style file (`gene<TAB>KO`, every row's gene id
  unique) was silently split into one fake "sample" per gene, because the
  `auto` heuristic only checked "more than one distinct first-column value" --
  which a per-gene file also satisfies. It now requires a sample id to
  actually *repeat* across rows before trusting the long-table reading; a
  file where every row's id is unique is read as a single sample instead
  (matching this tool's one-file-per-genome primary use case).
- New `--input-format long` forces long-table parsing for the rare genuine
  case where every sample happens to contribute exactly one row (structurally
  indistinguishable from the case above -- the caller must say so explicitly).
- Loading a directory where two files share a stem across extensions (e.g.
  `sample.txt` and `sample.tsv`) silently let the later one overwrite the
  earlier one's KOs with no warning. Now raises a clear error.

### Sample sheets
- A blank `sample_id` cell was silently accepted as a real (empty-string)
  sample. Now rejected with the offending row number(s).
- A blank `annotation_file` cell resolved, via `Path(base) / ""` being a
  no-op join, to the sheet's own base directory -- silently unioning every
  file in that directory into the row's KO set. Now rejected before it can
  load anything.
- `mpph validate` previously only checked that the `annotation_file` column
  existed, not that its values were non-blank or pointed at real files. It
  now reports both, so a bad sheet is caught before a real run instead of
  producing a quietly-wrong sample.

## 3.5.2

Correctness hotfix for ordination/PERMANOVA on degenerate data. Reproduced
with a failing case first, now has a regression test.

### Ordination
- `pcoa`'s `zero_distance_pairs` diagnostic counted `np.triu(d, 1) == 0`,
  which zeroes the diagonal *and* the whole lower triangle before comparing
  -- those artificial zeros were counted too, so 3 pairwise-*distinct*
  samples could report 6 "zero-distance pairs" that don't exist. Now counts
  ties only within the true strict upper triangle.
- `plot_ordination` crashed (`IndexError`) on a rank-one PCoA solution (e.g.
  exactly 2 samples, or other degenerate distance structures give only one
  positive eigenvalue). It now plots PCo1 with y fixed at 0 and labels the
  y-axis "PCo2 unavailable (rank-one solution)" instead of fabricating a
  second axis.

### PERMANOVA
- All-identical samples (every pairwise distance zero) produced `pseudo_F:
  nan` and a numerically meaningless but plausible-looking p-value (plus a
  wall of `RuntimeWarning: invalid value encountered in scalar divide`).
  Now raises `ValueError` before permuting, matching the precedent already
  set by `pcoa`'s own all-zero-distance guard.
- `permutations <= 0` silently returned a degenerate result (`permutations=0`
  gave p=1.0 with no test ever run; negative values gave a **negative**
  p-value). Now rejected with a clear `ValueError`.
- Non-finite distances (NaN/Inf, e.g. from a metric misapplied to the input)
  are now rejected before the permutation loop instead of propagating into
  the pseudo-F statistic silently.

## 3.5.1

Correctness hotfix for Newick tree export. Reproduced with a failing case
first, now has a regression test.

### Tree export
- `plot_matrix`'s returned linkage (`row_link`/`col_link`) indexes leaves in
  the *pre*-clustering row/column order, but the CLI's Newick export was
  passing it together with `row_labels`/`col_labels` -- the *post*-dendrogram,
  reordered display labels. The topology was correct but names could attach to
  the wrong leaves (e.g. two genuinely distant organisms shown as each
  other's closest relative, and vice versa). `plot_matrix` now also returns
  `row_link_labels`/`col_link_labels` (labels in the same order the linkage
  was computed on), and `organism_tree.nwk`/`feature_tree.nwk` export uses
  those instead.

## 3.5.0

New subcommand: **`mpph gsea`** — rank-based enrichment ("GSEAPreranked" style,
Subramanian et al. 2005), complementing `mpph enrich`'s hypergeometric ORA.
Every gene is scored and ranked (no threshold chosen to define a study set
first), and categories are tested for skewing toward either end of the ranking.

- `--ranked-list FILE`: bring your own `gene<TAB>score` ranking (e.g. a
  DESeq2/edgeR/limma statistic) -- the recommended, rigorous path.
- `--expression FILE --metadata FILE --group-column X --group-a A --group-b B`:
  rank genes automatically with `signal2noise` (the original GSEA paper's
  statistic) or `log2fc`. Explicitly documented as simple/transparent, **not a
  substitute for a dedicated DE tool**.
- `--ontology kegg-pathway` / `kegg-module` / `go`, same category-membership
  machinery as `enrich` (GO again needs a user-supplied gene-to-GO mapping).
- Weighted running-sum statistic (`--weight`, default 1.0 = standard GSEA
  weighting), significance by gene-set permutation batched per distinct
  category size (`--permutations`), NES, BH-FDR q-values, and leading-edge
  genes per category.
- Unlike `enrich`, **no background is required and no gene is dropped** --
  every ranked gene counts as a "miss" in the running sum (documented as a
  deliberate, important difference from ORA's annotated-universe restriction).
- Two new plots: `plot_gsea_running` (the classic running-enrichment-score
  plot with a hit rug and ranking-metric bars) and `plot_gsea_summary` (top
  hits by NES, coloured by enrichment direction and significance).
- New public API: `load_ranked_list`, `load_expression_matrix`,
  `rank_from_expression`, `enrichment_score`, `gsea_analysis`.
- 21 new tests: the running-sum statistic cross-checked on synthetic rankings
  with known top/bottom/scattered gene sets, leading-edge extraction, the
  expression-ranking statistics, and an offline CLI end-to-end test (KEGG
  pathway, KEGG module, GO, and expression-based ranking) served from
  committed cache fixtures.
- New example, computed live against real KEGG data: every KO across 25
  *Synechococcus* genomes ranked by (marine prevalence − freshwater
  prevalence), tested against KEGG modules. Top hit is "Sulfate-sulfur
  assimilation" (raw p≈0.001, correct ecological direction — marine water is
  sulfate-rich, freshwater is not) -- shown deliberately even though it does
  not clear q<0.05 across all 131 modules tested, as an honest example of
  reading GSEA output past a single threshold.
- docs/methods.md, docs/subcommands.md, docs/outputs.md, README updated.

## 3.4.0

New feature: **`mpph enrich`** — hypergeometric over-representation (ORA) of a
gene/KO study set against KEGG pathways, KEGG modules, or GO terms.

- `--ontology kegg-pathway` / `kegg-module`: category membership from the
  global (non-organism-specific) `link/pathway/ko` / `link/module/ko` KEGG
  endpoints, so results don't depend on which organism the study set came from.
  `--background-organism CODE` is a shortcut for "this organism's full KO
  complement" (built on the existing `organism_kos()`).
- `--ontology go`: GO enrichment via a **user-supplied** gene-to-GO mapping (a
  long table, or an eggNOG-mapper `.annotations` file's `GO_terms` column) --
  KEGG itself has no GO annotations, so this is bring-your-own-mapping, not a
  new external data source.
- One-sided hypergeometric test (over-representation only) with BH-FDR
  correction (reusing the existing `benjamini_hochberg`); both study and
  background are restricted to the annotated universe before testing, and
  dropped/used counts are reported.
- A ranked bar chart (`plot_enrichment`) of the top hits, coloured by
  significance.
- New public API: `fetch_ko_pathway_membership`, `fetch_ko_module_membership`,
  `fetch_pathway_names`, `hypergeometric_enrichment`, `invert_membership`,
  `load_gene_go_map`.
- 21 new tests: pure hypergeometric math (cross-checked directly against
  `scipy.stats.hypergeom`), GO/id-list file parsing, and an offline CLI
  end-to-end test served from committed KEGG cache fixtures. Also validated
  live against the real KEGG API (see the new example below).
- New example: KEGG pathway enrichment of the 190 genes unique to
  *Prochlorococcus* MIT 9313 (low-light ecotype) versus AS9601 (high-light,
  streamlined) — recovers "Photosynthesis - antenna proteins" as a top hit,
  the published low-light antenna-gene expansion in this genus.

## 3.3.1

Correctness and security hotfix. Every item below was reproduced with a failing
case first and now has a regression test.

### Missing data (NaN means *unknown*, never absent)
- `benjamini_hochberg`: a single NaN p-value no longer turns **every** q-value
  into NaN; only finite p-values are corrected and set the rank denominator.
- `differential_features`: NaN is dropped per feature and per group, so
  prevalences use a **known-only denominator**; reports
  `n_*_total/known/unknown` and a `warning`, and skips untestable features
  (`--min-known-per-group`) instead of silently treating unknown as absent.
- `pan_classify`: known-only prevalence + `n_present/n_known/n_unknown/
  known_fraction`; features below `min_known_fraction` are labelled
  `insufficient-data` rather than mis-classified as core/shell/cloud.

### Statistics
- `pcoa` now raises a clear error when no positive axis exists (all-zero
  distances) instead of returning an empty result that crashed `ordination`,
  and returns diagnostics (positive/negative eigenvalues, negative fraction).
- `permanova` excludes samples whose group label is missing (NaN was previously
  counted as its own group) and rejects groups with <2 samples; writes
  `*_permanova.csv`.
- `accumulation_curve` is documented as **permutation/rarefaction**, not
  bootstrap, and now reports SD and a 95% interval (plotted as a band).

### Security
- HTML report: the payload is embedded in a `<script type="application/json">`
  block with `<` escaped (no `</script>` breakout), NaN is serialised as null,
  and all user-controlled text is written via `textContent`/DOM construction
  instead of `innerHTML`. A CSP is set as defence in depth.

### Tree comparison
- `robinson_foulds` no longer loses a root-level singleton outgroup from the
  shared leaf set (identical `(a,(b,(c,d)))` reported 3 of 4 leaves); leaf sets
  are parsed from the tree and non-overlapping leaves are reported.

### Provenance
- Manifest records `prevalence_state`, `complete_threshold` and
  `score_semantics`, and uses the argv actually passed (correct for
  `main(argv=...)`).
- Examples regenerated at 3.3.1 with relative paths — the previous manifests
  leaked an absolute local build path.

## 3.3.0

- **Completeness states**: `--prevalence-state {any,complete}` +
  `--complete-threshold` — filter on how often a module is *fully complete*,
  not merely detectable.
- **More built-in trait panels**: `respiration` (terminal oxidases) and
  `carbon_fixation` (Calvin / rTCA / Wood-Ljungdahl / 3-HP), alongside
  `biogeochemistry`.
- **Analysis plots**: `compare` now writes a volcano plot; `pan` writes a
  pan-class prevalence histogram and a pan/core accumulation curve.
- **Trait figures** are labelled "trait completeness" / "trait category".
- **`docs/`**: methods, input formats, subcommands, and output schema.

## 3.2.0

Turns MPPH from a heatmap tool into a functional-profile analysis toolkit. The
CLI is now subcommand-based; `mpph <taxon> ...` still works as a shortcut for
`mpph run`.

### New analysis
- **Evidence-aware scoring** (`mpph.modules.evaluate_module`): per-step
  matched/missing KOs, complete/partial/absent/unknown state, unresolved
  references, parser status. Surfaced by **`mpph explain`**.
- **`mpph compare`** — differential features between two groups (Fisher exact /
  Mann-Whitney U, odds ratio / Cliff's delta, BH-FDR).
- **`mpph pan`** — core / soft-core / shell / cloud classes + bootstrap
  pan/core accumulation curve.
- **`mpph ordination`** — PCoA (classical MDS) with a scatter plot, and
  PERMANOVA when a grouping is supplied.
- **`mpph community`** — pairwise metabolic complementarity (modules completed
  by a union of organisms that neither completes alone).
- **`mpph traits`** — score JSON/YAML metabolic-trait panels; ships a built-in
  biogeochemistry marker panel (N/S/C/methane cycles, photosynthesis).
- **`mpph report`** — self-contained interactive HTML report (searchable colour
  heatmap, QC and manifest tabs).
- **`mpph validate`** — sample-sheet checks.

### New input / metadata
- `mpph.samplesheet`: reproducible sample-sheet input; CheckM2 and GTDB-Tk
  metadata import.
- `--input-format eggnog` reads the eggNOG-mapper `KEGG_ko` column.

### Public API
- The KEGG REST client and data-access functions are exported from the
  top-level package (see the README "Public Python API").

### Comparison / trees
- `mpph.treecompare`: Robinson-Foulds distance to a reference tree and
  feature-bootstrap clade support.

## 3.1.1
- Expose the KEGG connection interface as public API; derive the request
  User-Agent from the package version.

## 3.1.0

### Scientific correctness
- **Module scorer**: nested module references (`M#####`) now resolve recursively
  (cycle-guarded) instead of always scoring 0; `--` placeholder steps are
  excluded from numerator and denominator; completeness scores only `Pathway`
  modules by default (`--all-modules` to include signature/reaction modules).
- **Presence mode** now restricts to the `Metabolism` BRITE top-level category
  by default (`--top-category`, `--all-categories`); the tool no longer mixes
  in Human-Diseases / Genetic-Information-Processing pathways silently.
- Distance metric now defaults per data type — Jaccard for binary presence,
  Euclidean for continuous completeness — and Jaccard/Dice are rejected on
  completeness data (added `braycurtis`, `cosine`).

### Reproducibility & QC
- Per-organism fetch errors are recorded rather than fatal; a `*_qc.csv` reports
  each organism's annotated-feature count and status.
- Clustered order is exported: `*_row_order.csv`, `*_col_order.csv`,
  `*_ordered_matrix.csv`, and `*_organism_tree.nwk` / `*_feature_tree.nwk`.
- Enriched `*_manifest.json`: full command, Python/platform, KEGG release,
  before/after feature counts, excluded organisms, and all filter settings.

### Usability
- `--match {word,prefix,exact,regex}` for taxon name matching (`--exact` kept as
  a deprecated alias for `prefix`).
- Input sources (taxon / `--codes` / `--user`) are now mutually exclusive.
- Two-object clustering works; row and column dendrograms are controlled
  independently. Feature ids are shown on the x-axis for small matrices.
- Stable category colours (fixed hues for common KEGG categories).
- Honest framing: figures/labels call binary results "pathway-map association".

### Housekeeping
- Removed the unused `seaborn` dependency; `requires-python >= 3.9`.
- Legacy prototype notebooks moved to `legacy/`.
- Added `CITATION.cff`, `DATA_SOURCES.md`, this changelog.

## 3.0.0
- Repackaged the single script into the installable `mpph` package with a
  console entry point; added KEGG module-completeness mode, user MAG/KO input,
  Newick export, on-disk caching, a pytest suite, and GitHub Actions CI.
