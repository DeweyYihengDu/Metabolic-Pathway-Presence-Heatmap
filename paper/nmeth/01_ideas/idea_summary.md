# Idea ranking — what could make MPPH a methods contribution

Bounded by `../literature_matrix.md`. No preliminary results are claimed
below; every idea is written as a question with a stated way to kill it.

## Ranking

| Rank | Idea | Novelty | Feasibility | Risk | One-line |
|---|---|---|---|---|---|
| 1 | **003** — sub-threshold ground-truth audit | 4 | 5 | 2 | Cheapest decisive experiment; establishes (or kills) the premise the whole programme rests on |
| 2 | **001** — clade-held-out re-assessment | 4 | 5 | 3 | Highest-impact standalone result; directly testable; applies a proven critique to a field that has not had it |
| 3 | **002** — evidence-resolved calibrated posterior | 5 | 3 | 4 | The actual methods contribution, but only viable if 003 and 001 land |
| 4 | **005** — phylogenetically structured MAG simulator | 4 | 4 | 3 | Strong supporting infrastructure; not a paper alone |
| 5 | **004** — do published differences survive joint correction | 3 | 3 | 4 | Best as the application section, not the headline; adversarial framing risk |

## Recommended first move: idea_003

Not because it is the most exciting — 002 is — but because it is the only one
that can be run immediately on data already sitting on Hermes, needs no new
modelling, and **its outcome determines whether 002 is worth attempting at
all**. It answers one factual question against curated ground truth: how much
true signal does the field's binary threshold discard?

- If precision immediately below threshold is near zero → KOfam's cut is well
  placed, the "detection" channel of 002 is empty, and 002 shrinks to two
  channels. Cheap to learn, and worth knowing before building anything.
- If precision is substantial and rises with taxonomic distance from the
  profile's source sequences → the field is discarding recoverable evidence
  precisely where environmental genomics operates, and 002 has an empirical
  foundation rather than a plausible story.

**Main failure mode for idea_003:** KEGG's reference is itself
similarity-derived and incomplete, so a sub-threshold hit scored as a false
positive may be a real function KEGG simply has not recorded. This biases
measured precision **downward**. The result is therefore a conservative lower
bound on recoverable signal — which is the safe direction for the argument,
but it must be stated as a bound and not as a point estimate, and the write-up
must not quietly upgrade it later.

## Why the others are deprioritised now

- **001** is nearly as cheap and arguably more publishable, but it needs GTDB
  download and a matched re-implementation of someone else's model, i.e. days
  not hours, and it carries a fairness obligation (run their released model,
  not my reconstruction of it). It should start as soon as 003 reports.
  Its mandatory control — size-matched random splits at each block size — is
  what separates a real finding from "smaller training set performs worse".
- **002** is the genuine contribution but has the highest chance of failing on
  identifiability. It should not be started before 003 and 001 supply its two
  premises.
- **005** is attractive and self-contained, but a simulation-only paper will
  not clear a top methods venue, and building a benchmark one's own tool wins
  on is a trap. It follows 002.
- **004** re-analyses other people's published conclusions. Scientifically
  fair, socially costly, and dependent on 002's machinery. It belongs in the
  final paper as the demonstration that the method changes conclusions, using
  the author's own PVC data as the primary case.

## Honest position on the target venue

The user asked for Nature Methods. Stated plainly, and not to be softened
later: MPPH as it stands today is a well-engineered, well-validated
integration toolkit whose two headline validations are *reproductions* of
KofamScan and clusterProfiler, and whose phylogenetic correction is a port of
an idea with four existing implementations (Scoary, treeWAS, pyseer, hogwash).
That is a strong Bioinformatics Applications Note. It is not a Nature Methods
paper, and submitting it as one would waste months.

The route that could reach a top methods venue is 003 → 001 → 002: establish
empirically that the field's binary calls discard real signal, show that its
accuracy claims are inflated by phylogenetic pseudoreplication, then offer a
calibrated alternative validated the way no incumbent has been. Even executed
well, Genome Biology / Nature Communications / Microbiome are the more likely
landing zones; Nature Methods is possible but should be treated as the
stretch outcome, not the plan.
