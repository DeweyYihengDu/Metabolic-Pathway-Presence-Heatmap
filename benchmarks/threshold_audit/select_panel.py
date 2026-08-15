"""Select a genome panel that spans a *range* of relatedness.

The neighbour-calibration question -- does fitting the confidence curve on a
closer relative calibrate a query better? -- cannot be answered by the seven
benchmark genomes, which are mutually maximally distant. It needs pairs at
every taxonomic rank: two strains of one species, two species of one genus,
two genera of one family, and so on up to different domains.

This picks such a panel from KEGG's own organism taxonomy (BRITE br08601 --
11,949 organisms with full lineage), which is one request rather than one per
organism, and requires only that an organism have a curated KO set to serve as
ground truth.

Selection is deliberately *structured*, not random: random sampling of KEGG
would return an overwhelmingly Proteobacterial panel with almost no
close-relative pairs at all, which is precisely the design flaw this whole
line of work is about. Instead, genera are chosen to contribute several
species each, so the close-pair cells are populated by construction.

    python select_panel.py --per-genus 4 --n-genera 20 --out panel.tsv
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from mpph.kegg import DEFAULT_CACHE, kegg_get, make_session

# br08601 lines look like:
#   AProkaryotes (10621)
#   B  Bacteria (10310)
#   C    Gammaproteobacteria (3122)
#   D      Escherichia (147)
#   E        eco  Escherichia coli K-12 MG1655
# Depth letter is positional rather than named, so the rank mapping is
# recorded here rather than guessed per line. Counts in parentheses are part
# of the label and are stripped.
_ENTRY = re.compile(r"^([A-G])\s*(.*)$")
_COUNT = re.compile(r"\s*\(\d+\)\s*$")


def parse_br08601(text: str) -> list[dict]:
    """Return one record per organism with its enclosing lineage."""
    lineage: dict[str, str] = {}
    out = []
    for raw in text.splitlines():
        m = _ENTRY.match(raw.rstrip())
        if not m:
            continue
        depth, rest = m.group(1), m.group(2).strip()
        rest = _COUNT.sub("", re.sub(r"</?b>", "", rest)).strip()
        if not rest:
            continue
        parts = rest.split(None, 1)
        # An organism line starts with a 3-4 letter KEGG org code followed by
        # a name; a rank line is just a name.
        if len(parts) == 2 and re.fullmatch(r"[a-z]{3,4}", parts[0]):
            out.append({"org": parts[0], "name": parts[1], **lineage})
        else:
            lineage = {k: v for k, v in lineage.items() if k < depth}
            lineage[depth] = rest
    return out


def genus_of(name: str) -> str:
    """First word of the binomial, skipping KEGG's placeholder prefixes.

    Without this, "Candidatus", "Endosymbiont" and "uncultured" become the
    largest "genera" in the panel and every within-genus pair they contribute
    is spurious -- organisms sharing only a nomenclatural placeholder are not
    close relatives.
    """
    if not name:
        return ""
    words = name.replace("[", "").replace("]", "").split()
    skip = {"Candidatus", "candidatus", "uncultured", "Uncultured",
            "Endosymbiont", "endosymbiont", "Candidate"}
    for w in words:
        if w not in skip:
            # A real genus is capitalised and alphabetic; "sp." / "bacterium"
            # style placeholders are not usable as a genus key.
            return w if w[:1].isupper() and w.isalpha() else ""
    return ""


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--per-genus", type=int, default=4,
                   help="organisms to take from each chosen genus")
    p.add_argument("--n-genera", type=int, default=20)
    p.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    p.add_argument("--out", required=True)
    args = p.parse_args()

    session = make_session()
    text = kegg_get(session, "get/br:br08601", Path(args.cache_dir))
    records = parse_br08601(text)
    print(f"{len(records)} organism(s) in br08601")

    # Key on (kingdom, genus), never genus alone. Genus names are not unique
    # across kingdoms -- KEGG lists both a bacterial *Bacillus* and an insect
    # *Bacillus* -- and pooling them would put a stick insect and a soil
    # bacterium in the same "within-genus" cell, which is exactly the kind of
    # pair this panel exists to measure correctly.
    by_genus: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in records:
        if not r.get("A", "").startswith("Prokaryotes"):
            continue      # this is the KOfam calibration question; eukaryotes
                          # are the separate tier already covered by sce/ath
        g = genus_of(r["name"])
        if g:
            by_genus[(r.get("A", ""), g)].append(r)

    # Genera with enough members to give within-genus pairs, ordered by size
    # so the panel is populated where the close-pair cells need it.
    candidates = sorted((k for k, v in by_genus.items()
                         if len(v) >= args.per_genus),
                        key=lambda k: -len(by_genus[k]))

    # Spread across high-level lineages rather than taking the N biggest
    # genera, which would be almost all Proteobacteria.
    chosen, seen_lineage = [], defaultdict(int)
    for key in candidates:
        top = by_genus[key][0].get("C", "?")
        if seen_lineage[top] >= 3:
            continue
        seen_lineage[top] += 1
        chosen.append(key)
        if len(chosen) >= args.n_genera:
            break

    panel = []
    for key in chosen:
        panel.extend(by_genus[key][:args.per_genus])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("org\tname\tkingdom\tdomain\tclade\tgenus\n")
        for r in panel:
            fh.write("\t".join([r["org"], r["name"], r.get("A", ""),
                                r.get("B", ""), r.get("C", ""),
                                genus_of(r["name"])]) + "\n")
    print(f"{len(panel)} organism(s) from {len(chosen)} genera across "
          f"{len(seen_lineage)} clades. Wrote {out}")
    for key in chosen:
        print(f"  {key[1]:22s} ({by_genus[key][0].get('C','?')})  "
              f"n={min(args.per_genus, len(by_genus[key]))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
