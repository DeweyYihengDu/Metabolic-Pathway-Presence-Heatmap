"""Fetch a KO presence matrix + taxonomy for the phylogenetic-leakage test.

The question (idea_001) is whether reported accuracy for genome-based
prediction of KEGG module presence survives evaluation on phylogenetically
held-out clades, or is inflated by pseudoreplication between train and test
sets drawn from a clade-redundant database.

Answering it needs KO vectors and taxonomy for a few thousand genomes -- and
nothing else. No proteomes, no HMM search, no GTDB: KEGG publishes the curated
KO set per organism (`link/ko/<org>`) and the full lineage (BRITE br08601)
directly. That makes this a network-bound job that can run alongside a
CPU-bound one.

**Politeness.** KEGG REST is a free public service. Requests are issued
serially with a delay, responses are cached on disk by `mpph.kegg` so a rerun
costs nothing, and the organism count is a parameter rather than "all of
KEGG". Interrupting and resuming is safe.

    python fetch_ko_matrix.py --n-organisms 3000 --out-dir ko_matrix/
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "threshold_audit"))

from select_panel import genus_of, parse_br08601

from mpph.kegg import DEFAULT_CACHE, kegg_get, make_session


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-organisms", type=int, default=3000)
    p.add_argument("--delay", type=float, default=0.34,
                   help="seconds between requests; KEGG asks for <=3/sec")
    p.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()

    session = make_session()
    cache = Path(args.cache_dir)
    out = Path(args.out_dir)
    (out / "ko").mkdir(parents=True, exist_ok=True)

    records = parse_br08601(kegg_get(session, "get/br:br08601", cache))
    prok = [r for r in records if r.get("A", "").startswith("Prokaryotes")]
    print(f"{len(records)} organisms in br08601, {len(prok)} prokaryotic")

    # Take every organism in taxonomic order rather than a random sample: the
    # point is to reproduce the clade redundancy a real training set has, and
    # random sampling would thin exactly the dense clades that cause the
    # problem being measured.
    chosen = prok[:args.n_organisms]

    tax_path = out / "taxonomy.tsv"
    with open(tax_path, "w", encoding="utf-8") as fh:
        fh.write("org\tname\tkingdom\tdomain\tclade\tgenus\n")
        fh.writelines("\t".join([r["org"], r["name"], r.get("A", ""),
                                r.get("B", ""), r.get("C", ""),
                                genus_of(r["name"])]) + "\n" for r in chosen)
    print(f"wrote {tax_path}")

    n_ok = n_skip = n_cached = 0
    for i, r in enumerate(chosen, 1):
        dest = out / "ko" / f"{r['org']}.tsv"
        if dest.exists():
            n_cached += 1
            continue
        try:
            text = kegg_get(session, f"link/ko/{r['org']}", cache)
        except Exception as exc:                      # noqa: BLE001
            print(f"\n  {r['org']}: {exc}")
            n_skip += 1
            continue
        if not text.strip():
            n_skip += 1
            dest.write_text("", encoding="utf-8")     # remember the miss
            continue
        dest.write_text(text, encoding="utf-8")
        n_ok += 1
        time.sleep(args.delay)
        if i % 25 == 0:
            print(f"\r  {i}/{len(chosen)} ok={n_ok} cached={n_cached} "
                  f"empty={n_skip}", end="", flush=True)
    print(f"\n{n_ok} fetched, {n_cached} already present, {n_skip} with no KO set")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
