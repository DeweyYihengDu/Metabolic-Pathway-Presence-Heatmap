"""Head-to-head comparison: mpph annotate vs. KofamScan vs. eggNOG-mapper.

Compares per-gene KO assignments from all three tools against each other and
against each organism's real, KEGG-curated ground truth (`link/ko/<org>`).
Not part of the installed `mpph` package -- this is benchmark tooling that
reads each tool's own output format and reports agreement/accuracy tables.

Usage:
    python compare_annotation.py --org eco \
        --mpph results/mpph_eco.tsv \
        --kofamscan results/kofamscan_eco.tsv \
        --emapper results/emapper_eco.emapper.annotations \
        --ground-truth ground_truth/link_ko_eco.tsv \
        --protein-id-map ground_truth/protein_id_map_eco.tsv \
        --n-genes 4300 \
        --out report/annotation_eco.csv
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

_KO = re.compile(r"K\d{5}")


def _strip_version(accession: str) -> str:
    """`NP_414542.1` -> `NP_414542`. The genome FASTA (and therefore every
    tool's own gene-id output) carries the RefSeq version suffix; KEGG's
    `conv` tables do not -- stripped here so both sides match."""
    return accession.rsplit(".", 1)[0] if re.search(r"\.\d+$", accession) else accession


def load_mapper_tsv(path: str | Path) -> dict[str, set[str]]:
    """`gene<TAB>K#####` (no header), the shape both `mpph annotate` and
    KofamScan `--format mapper` write. Returns {gene_id: {KO, ...}}.

    KofamScan's own mapper format writes one line per *query* gene, even
    when it found no significant hit -- such a line is just the bare gene
    id with no tab/KO at all (confirmed against real KofamScan output: ~23%
    of E. coli's lines are exactly this). `mpph annotate` instead omits a
    gene entirely when it has no assignment. Both conventions mean "no KO
    call for this gene" and are treated identically here -- skipped, not
    an error.
    """
    pairs: dict[str, set[str]] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            continue  # gene with no significant hit (KofamScan's convention)
        gene, ko = fields[:2]
        pairs.setdefault(_strip_version(gene), set()).add(ko)
    return pairs


def load_emapper_annotations(path: str | Path) -> dict[str, set[str]]:
    """eggNOG-mapper's `.annotations` file: `#query` id column + `KEGG_ko`
    column (comma-separated `ko:K#####` values, or `-` for none)."""
    header: list[str] | None = None
    idx_query = idx_ko = None
    pairs: dict[str, set[str]] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("##") or not line.strip():
            continue
        fields = line.lstrip("#").split("\t")
        if header is None:
            header = fields
            idx_query = 0  # first column is the query/gene id in every emapper version
            idx_ko = header.index("KEGG_ko")
            continue
        if line.startswith("#"):
            continue
        gene = fields[idx_query]
        kos = set(_KO.findall(fields[idx_ko])) if idx_ko < len(fields) else set()
        if kos:
            pairs[_strip_version(gene)] = kos
    return pairs


def load_protein_id_map(path: str | Path, org: str) -> dict[str, str]:
    """KEGG's `conv/<org>/ncbi-proteinid` -- `ncbi-proteinid:NP_414542<TAB>
    org:b0001` -- returns {KEGG_gene_id: ncbi_protein_id} (b0001 -> NP_414542,
    no version suffix). Needed because KEGG's `link/ko/<org>` ground truth is
    keyed by KEGG's own gene id, but every annotation tool here is keyed by
    whatever the input FASTA's headers were -- NCBI RefSeq protein accessions
    -- a genuinely different id namespace for the same genes, not a
    tool-accuracy difference.
    """
    prefix = f"{org}:"
    mapping: dict[str, str] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        protein_field, gene_field = line.split("\t")[:2]
        protein_id = protein_field.split(":", 1)[-1]
        gene_id = gene_field.removeprefix(prefix)
        mapping[gene_id] = protein_id
    return mapping


def load_ground_truth(path: str | Path, org: str,
                      protein_id_map: dict[str, str] | None = None) -> dict[str, set[str]]:
    """KEGG's own `link/ko/<org>` -- `org:gene<TAB>ko:K#####`.

    If `protein_id_map` is given, KEGG gene ids are translated to NCBI
    protein accessions (see `load_protein_id_map`) so this lines up with the
    id namespace every annotation tool's own output actually uses. A KEGG
    gene with no corresponding NCBI protein id in the map (rare, but real --
    e.g. RNA genes with no protein product) is dropped rather than silently
    kept under its untranslatable KEGG id, which would never match anything.
    """
    pairs: dict[str, set[str]] = {}
    prefix = f"{org}:"
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        gene_field, ko_field = line.split("\t")[:2]
        gene = gene_field.removeprefix(prefix)
        ko = ko_field.split(":", 1)[-1]
        if protein_id_map is not None:
            gene = protein_id_map.get(gene)
            if gene is None:
                continue
        pairs.setdefault(gene, set()).add(ko)
    return pairs


def _pairs(gene_to_kos: dict[str, set[str]]) -> set[tuple[str, str]]:
    """Flatten {gene: {KO, ...}} to a set of (gene, KO) pairs -- the natural
    unit for precision/recall when a gene can have more than one true KO."""
    return {(gene, ko) for gene, kos in gene_to_kos.items() for ko in kos}


def precision_recall_f1(predicted: dict[str, set[str]],
                        truth: dict[str, set[str]]) -> dict[str, float]:
    pred_pairs, true_pairs = _pairs(predicted), _pairs(truth)
    tp = len(pred_pairs & true_pairs)
    fp = len(pred_pairs - true_pairs)
    fn = len(true_pairs - pred_pairs)
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    f1 = (2 * precision * recall / (precision + recall)
          if precision and recall and (precision + recall) else float("nan"))
    return {"tp": tp, "fp": fp, "fn": fn,
            "precision": precision, "recall": recall, "f1": f1}


def pairwise_agreement(a: dict[str, set[str]], b: dict[str, set[str]]) -> dict[str, float]:
    """Jaccard overlap of (gene, KO) pairs between two tools -- symmetric,
    unlike precision/recall which need a designated ground truth."""
    pa, pb = _pairs(a), _pairs(b)
    union = pa | pb
    jaccard = len(pa & pb) / len(union) if union else float("nan")
    return {"n_a": len(pa), "n_b": len(pb), "n_shared": len(pa & pb), "jaccard": jaccard}


def coverage(gene_to_kos: dict[str, set[str]], n_genes: int) -> float:
    """Fraction of the genome's genes this tool assigned >=1 KO to."""
    return len(gene_to_kos) / n_genes if n_genes else float("nan")


def build_report(org: str, n_genes: int,
                 mpph: dict, kofamscan: dict, emapper: dict,
                 truth: dict | None) -> pd.DataFrame:
    """`truth=None` (a genome with no KEGG-curated KO assignments -- e.g. a
    MAG assembled from a metagenome) reports coverage and pairwise tool
    agreement only. Precision/recall are simply undefined there: absence of
    a reference is not evidence a call is wrong, and scoring against a
    non-existent truth set would manufacture numbers rather than measure
    anything."""
    tools = {"mpph_annotate": mpph, "kofamscan": kofamscan, "eggnog_mapper": emapper}
    rows = []
    for name, calls in tools.items():
        row = {"organism": org, "tool": name, "coverage": coverage(calls, n_genes),
               "n_calls": len(_pairs(calls))}
        if truth is not None:
            row.update(precision_recall_f1(calls, truth))
        rows.append(row)
    accuracy = pd.DataFrame(rows)

    pair_rows = []
    names = list(tools)
    for i, name_a in enumerate(names):
        for name_b in names[i + 1:]:
            agreement = pairwise_agreement(tools[name_a], tools[name_b])
            pair_rows.append({"organism": org, "tool_a": name_a, "tool_b": name_b,
                              **agreement})
    agreement_df = pd.DataFrame(pair_rows)
    return accuracy, agreement_df


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--org", required=True)
    p.add_argument("--n-genes", type=int, required=True)
    p.add_argument("--mpph", required=True)
    p.add_argument("--kofamscan", required=True)
    p.add_argument("--emapper", required=True)
    p.add_argument("--ground-truth", default=None,
                   help="KEGG link/ko/<org> table. Omit for a genome with no "
                        "KEGG-curated annotation (e.g. a MAG) -- coverage and "
                        "pairwise tool agreement are still reported, "
                        "precision/recall are not (there is nothing to score "
                        "against).")
    p.add_argument("--protein-id-map", default=None,
                   help="KEGG conv/<org>/ncbi-proteinid table -- translates "
                        "the ground truth's KEGG gene ids to the NCBI "
                        "protein accessions the input FASTA (and therefore "
                        "every tool's own output) actually uses. Omit only "
                        "if your genome FASTA already uses KEGG gene ids.")
    p.add_argument("--out", required=True, help="Accuracy table CSV path.")
    p.add_argument("--out-agreement", default=None,
                   help="Pairwise agreement table CSV (default: <out>_agreement.csv)")
    args = p.parse_args()

    mpph = load_mapper_tsv(args.mpph)
    kofamscan = load_mapper_tsv(args.kofamscan)
    emapper = load_emapper_annotations(args.emapper)
    protein_id_map = (load_protein_id_map(args.protein_id_map, args.org)
                      if args.protein_id_map else None)
    truth = (load_ground_truth(args.ground_truth, args.org, protein_id_map)
             if args.ground_truth else None)

    accuracy, agreement = build_report(args.org, args.n_genes, mpph, kofamscan, emapper, truth)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    accuracy.to_csv(out, index=False)
    agreement_path = Path(args.out_agreement or out.with_name(f"{out.stem}_agreement.csv"))
    agreement.to_csv(agreement_path, index=False)

    print(accuracy.to_string(index=False))
    print()
    print(agreement.to_string(index=False))
    print(f"\nWrote {out} and {agreement_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
