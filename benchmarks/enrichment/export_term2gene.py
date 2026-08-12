"""Export mpph's own KEGG pathway <-> KO membership snapshot as a
clusterProfiler-ready TERM2GENE table (pathway_id, KO_id columns).

This is the mechanism that makes the mpph-vs-clusterProfiler comparison fair:
both tools are tested against the *identical* category definitions (the same
KEGG release, fetched once here), so any difference in results reflects the
statistical implementation, not a different underlying database snapshot.
Reuses mpph's own fetch/filter functions directly -- no reimplementation.

Usage:
    python export_term2gene.py --top-category Metabolism --out term2gene.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # repo root, for `import mpph`

from mpph.cli import _restrict_to_top_category
from mpph.enrichment import (
    fetch_ko_pathway_membership,
    fetch_pathway_names,
    invert_membership,
)
from mpph.kegg import DEFAULT_CACHE, make_session


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--top-category", default="Metabolism",
                   help="Restrict to this BRITE top-level category, or omit "
                        "to export every KEGG pathway.")
    p.add_argument("--all-categories", action="store_true")
    p.add_argument("--out", required=True)
    p.add_argument("--out-names", default=None,
                   help="Optional pathway_id<TAB>name table (TERM2NAME).")
    p.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    args = p.parse_args()

    session = make_session()
    cache_dir = Path(args.cache_dir)

    item_to_cats = fetch_ko_pathway_membership(session, cache_dir)
    category_to_items = invert_membership(item_to_cats)
    top_category = None if args.all_categories else args.top_category
    category_to_items, _ = _restrict_to_top_category(
        category_to_items, session, cache_dir, refresh=False,
        top_category=top_category)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        for pathway_id, kos in sorted(category_to_items.items()):
            fh.writelines(f"{pathway_id}\t{ko}\n" for ko in sorted(kos))
    n_pairs = sum(len(v) for v in category_to_items.values())
    print(f"{len(category_to_items)} pathway(s), {n_pairs} (pathway, KO) pairs. Wrote {out}")

    if args.out_names:
        names = fetch_pathway_names(session, cache_dir)
        names_out = Path(args.out_names)
        with open(names_out, "w", encoding="utf-8") as fh:
            for pathway_id in sorted(category_to_items):
                fh.write(f"{pathway_id}\t{names.get(pathway_id, pathway_id)}\n")
        print(f"Wrote {names_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
