"""Retain KOfam hits that fall *below* each KO's adaptive threshold.

Every tool in this space -- KofamScan, mpph annotate, DRAM, MicrobeAnnotator,
anvi'o -- converts the HMM search into a binary call at that threshold and
discards everything under it. This script keeps the discarded band so its
precision against KEGG's own curated assignments can be measured
(`score_precision.py`).

The question being answered is factual, not modelled: among (gene, KO) pairs
scoring just below threshold, what fraction does KEGG itself assign?

Scoring mirrors `mpph.kofam.is_significant_hit` exactly -- best-domain score
for `score_type == "domain"` KOs, full-sequence score otherwise -- so the
delta computed here is the same quantity the real significance test uses,
just not thresholded. A hit with delta >= 0 is precisely what mpph would have
called; delta < 0 is the discarded band.

Only hits with `delta >= -DELTA_FLOOR` are written. At T=0 the search returns
an enormous number of near-zero-score hits that are of no interest and would
make the output unmanageable; the floor is generous enough to see where
precision decays to the background rate.

    python collect_subthreshold.py --fasta genomes/eco.faa \
        --kofam-db databases/kofam --cpus 24 --out subthreshold_eco.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # repo root

from mpph.kofam import discover_profiles, load_ko_thresholds

DELTA_FLOOR = 60.0


def collect(fasta_path: Path, kofam_db_dir: Path, *, cpus: int,
            delta_floor: float = DELTA_FLOOR):
    """Yield (gene_id, ko_id, score_used, threshold, delta, score_type)."""
    import pyhmmer

    entries = load_ko_thresholds(kofam_db_dir)
    profile_paths = discover_profiles(kofam_db_dir)
    alphabet = pyhmmer.easel.Alphabet.amino()

    def _iter_profile_hmms():
        for p in profile_paths:
            with pyhmmer.plan7.HMMFile(p) as hmm_file:
                yield from hmm_file

    def _name(x):
        return x.decode() if isinstance(x, bytes) else x

    with pyhmmer.easel.SequenceFile(fasta_path, digital=True,
                                    alphabet=alphabet) as sf:
        sequences = sf.read_block()

    for hits in pyhmmer.hmmsearch(_iter_profile_hmms(), sequences,
                                  cpus=cpus, T=0):
        ko = _name(hits.query.name)
        entry = entries.get(ko)
        if entry is None or not entry.assignable:
            continue  # no threshold in ko_list -- can never be assigned at all
        for hit in hits:
            if entry.score_type == "domain":
                if not len(hit.domains):
                    continue
                score_used = hit.best_domain.score
            else:
                score_used = hit.score
            delta = score_used - entry.threshold
            if delta < -delta_floor:
                continue
            yield (_name(hit.name), ko, score_used, entry.threshold, delta,
                   entry.score_type)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta", required=True)
    p.add_argument("--kofam-db", required=True)
    p.add_argument("--cpus", type=int, default=0)
    p.add_argument("--delta-floor", type=float, default=DELTA_FLOOR)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("gene\tko\tscore\tthreshold\tdelta\tscore_type\n")
        for gene, ko, score, thr, delta, stype in collect(
                Path(args.fasta), Path(args.kofam_db), cpus=args.cpus,
                delta_floor=args.delta_floor):
            fh.write(f"{gene}\t{ko}\t{score:.4f}\t{thr:.4f}\t{delta:.4f}\t{stype}\n")
            n += 1
    print(f"{n} hit(s) with delta >= -{args.delta_floor}. Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
