"""Command-line interface for MPPH."""
from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

from . import __version__
from .kegg import DEFAULT_CACHE, kegg_release, make_session
from .matrix import (
    build_completeness_matrix,
    build_presence_matrix,
    filter_matrix,
)
from .modules import fetch_module_definitions, list_modules, organism_kos
from .organisms import (
    list_genomes,
    read_code_file,
    select_by_codes,
    select_by_taxon,
)
from .pathways import OVERVIEW_CATEGORY, fetch_pathway_categories, get_pathways
from .plot import plot_matrix
from .tree import linkage_to_newick

DEFAULT_CACHE_STR = str(DEFAULT_CACHE)


def run(args: argparse.Namespace) -> int:
    cache_dir = None if args.no_cache else Path(args.cache_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    session = make_session()

    from_user = args.user is not None
    mode = "completeness" if (args.completeness or from_user) else "presence"
    # Distance defaults suit the data: Jaccard for binary presence, Euclidean
    # for continuous completeness. Reject metrics that would binarise it.
    metric = args.metric or ("euclidean" if mode == "completeness" else "jaccard")
    if mode == "completeness" and metric in {"jaccard", "dice"}:
        print(f"Error: metric {metric!r} binarises continuous completeness; "
              "use euclidean, braycurtis, or cosine.", file=sys.stderr)
        return 2
    match = "prefix" if args.exact else args.match
    label_src = args.user if from_user else (args.taxon or args.codes or "taxon")
    slug = re.sub(r"[^0-9A-Za-z]+", "_", Path(str(label_src)).stem).strip("_") or "mpph"
    failed: list[dict] = []

    # --- gather each organism's feature content ------------------------------
    if from_user:
        from .userdata import load_user_kos
        print(f"[1/4] Loading user KO annotations from {args.user} ...", flush=True)
        org_kos = load_user_kos(args.user)
        organisms = [(name, name) for name in org_kos]
        print(f"      loaded {len(org_kos)} sample(s)")
    else:
        print("[1/4] Fetching KEGG genome list ...", flush=True)
        genomes = list_genomes(session, cache_dir, refresh=args.refresh)
        if args.codes:
            organisms = select_by_codes(genomes, read_code_file(args.codes))
        else:
            organisms = select_by_taxon(genomes, args.taxon, match)
        if not organisms:
            print(f"No KEGG genomes matched {label_src!r}.", file=sys.stderr)
            print("Try a broader name, a different --match, or check spelling.",
                  file=sys.stderr)
            return 2
        print(f"      matched {len(organisms)} genome(s)")

    # --- build the matrix (per-organism errors are recorded, not fatal) ------
    top_categories: dict[str, str] = {}
    if mode == "completeness":
        print("[2/4] Fetching module definitions + organism KOs ...", flush=True)
        modules = list_modules(session, cache_dir, refresh=args.refresh)
        module_defs = fetch_module_definitions(
            session, sorted(modules), cache_dir, refresh=args.refresh,
            progress=lambda n, t: print(f"      module defs {n}/{t}",
                                        end="\r", flush=True),
        )
        print()
        if not args.all_modules:
            module_defs = {m: v for m, v in module_defs.items() if v[2] == "Pathway"}
            print(f"      {len(module_defs)} pathway modules "
                  "(use --all-modules to include signature/reaction modules)")
        categories = {mid: cat for mid, (_d, cat, _t) in module_defs.items()}
        feature_names = modules
        if from_user:
            kos = org_kos
        else:
            kos = {}
            for i, (code, name) in enumerate(organisms, 1):
                try:
                    kos[f"{name} ({code})"] = organism_kos(
                        session, code, cache_dir, refresh=args.refresh)
                except requests.RequestException as exc:
                    failed.append({"organism": f"{name} ({code})",
                                   "reason": f"ko_fetch_error: {exc}"})
                print(f"      KOs {i}/{len(organisms)}", end="\r", flush=True)
            print()
        df, _ = build_completeness_matrix(kos, module_defs)
    else:
        print(f"[2/4] Fetching pathways for {len(organisms)} genome(s) ...",
              flush=True)
        org_pathways = {}
        for i, (code, name) in enumerate(organisms, 1):
            try:
                org_pathways[f"{name} ({code})"] = get_pathways(
                    session, code, cache_dir, refresh=args.refresh)
            except requests.RequestException as exc:
                failed.append({"organism": f"{name} ({code})",
                               "reason": f"pathway_fetch_error: {exc}"})
                continue
            print(f"      [{i}/{len(organisms)}] {name} ({code})", flush=True)
        cat_full = fetch_pathway_categories(session, cache_dir, refresh=args.refresh)
        categories = {m: b for m, (_a, b) in cat_full.items()}
        top_categories = {m: a for m, (a, _b) in cat_full.items()}
        df, feature_names = build_presence_matrix(org_pathways)

    if failed:
        print(f"      warning: {len(failed)} organism(s) failed to fetch "
              "(recorded in the QC table)", file=sys.stderr)
    annotated_counts = (df > 0).sum(axis=1).to_dict()  # for the QC report

    # --- filter --------------------------------------------------------------
    # Order: restrict to a top-level category, drop aggregate overview maps,
    # then organisms left with no real features, then prevalence/core filtering
    # -- one all-absent genome would otherwise cap every prevalence below 1.0.
    print("[3/4] Filtering matrix ...", flush=True)
    n_features_raw = int(df.shape[1])
    if mode == "presence" and not args.all_categories:
        keep = {m for m, a in top_categories.items() if a == args.top_category}
        df = df.loc[:, [c for c in df.columns if c in keep]]
    if mode == "presence" and not args.keep_overview:
        overview = {m for m, c in categories.items() if c == OVERVIEW_CATEGORY}
        df = df.drop(columns=[c for c in df.columns if c in overview])
    excluded_empty: list[str] = []
    if not args.keep_empty:
        excluded_empty = list(df.index[df.sum(axis=1) == 0])
        if excluded_empty:
            print(f"      dropping {len(excluded_empty)} organism(s) with no "
                  f"retained features: {', '.join(excluded_empty)}", flush=True)
            df = df.drop(index=excluded_empty)
        if df.shape[0] < 1:
            print("All organisms were empty after filtering. Nothing to plot.",
                  file=sys.stderr)
            return 3
    df = filter_matrix(df, args.min_prevalence, args.max_prevalence,
                       args.drop_core, present_threshold=args.present_threshold)
    if df.shape[1] == 0:
        print("No features left after filtering. Loosen the filters.",
              file=sys.stderr)
        return 3
    print(f"      matrix: {df.shape[0]} organisms x {df.shape[1]} features")

    # --- write tables --------------------------------------------------------
    matrix_csv = outdir / f"{slug}_matrix.csv"
    features_csv = outdir / f"{slug}_features.csv"
    qc_csv = outdir / f"{slug}_qc.csv"
    df.to_csv(matrix_csv, index_label="organism")
    kept = list(df.columns)
    pd.DataFrame({
        "feature_id": kept,
        "name": [feature_names.get(k, k) if feature_names else k for k in kept],
        "category": [categories.get(k, "Other") for k in kept],
        "top_category": [top_categories.get(k, "") for k in kept],
    }).to_csv(features_csv, index=False)

    qc_rows = [{"organism": o, "n_annotated_features": int(n),
                "status": "excluded_empty" if o in excluded_empty else "included"}
               for o, n in annotated_counts.items()]
    qc_rows += [{"organism": f["organism"], "n_annotated_features": 0,
                 "status": "fetch_failed"} for f in failed]
    pd.DataFrame(qc_rows).to_csv(qc_csv, index=False)

    # --- figure(s) -----------------------------------------------------------
    print("[4/4] Rendering heatmap ...", flush=True)
    verb = "pathway-map association" if mode == "presence" else "module completeness"
    title = f"KEGG {verb} · {label_src}"
    subtitle = (f"{df.shape[0]} organisms × {df.shape[1]} "
                + ("KEGG modules" if mode == "completeness" else "KEGG pathways")
                + (f"  ·  UPGMA ({metric})" if args.cluster else ""))
    figures, layout = [], {}
    for fmt in args.format:
        fig_path = outdir / f"{slug}_heatmap.{fmt}"
        layout = plot_matrix(df, categories, fig_path, title, subtitle,
                             cluster=args.cluster, mode=mode, metric=metric)
        figures.append(fig_path.name)

    # --- clustered order + trees (traceability) ------------------------------
    ordered_outputs: dict[str, str] = {}
    if layout.get("row_labels"):
        ordered = df.loc[layout["row_labels"], layout["col_labels"]]
        f_ord = outdir / f"{slug}_ordered_matrix.csv"
        ordered.to_csv(f_ord, index_label="organism")
        ro = outdir / f"{slug}_row_order.csv"
        pd.DataFrame({"position": range(1, len(layout["row_labels"]) + 1),
                      "organism": layout["row_labels"]}).to_csv(ro, index=False)
        co = outdir / f"{slug}_col_order.csv"
        pd.DataFrame({
            "position": range(1, len(layout["col_labels"]) + 1),
            "feature_id": layout["col_labels"],
            "name": [feature_names.get(k, k) if feature_names else k
                     for k in layout["col_labels"]],
        }).to_csv(co, index=False)
        ordered_outputs = {"ordered_matrix": f_ord.name, "row_order": ro.name,
                           "col_order": co.name}

    # --- Newick trees --------------------------------------------------------
    tree_files: dict[str, str] = {}
    if args.newick:
        if layout.get("row_link") is not None:
            f = outdir / f"{slug}_organism_tree.nwk"
            f.write_text(linkage_to_newick(layout["row_link"],
                                           layout["row_labels"]), encoding="utf-8")
            tree_files["organism_tree"] = f.name
        if layout.get("col_link") is not None:
            f = outdir / f"{slug}_feature_tree.nwk"
            f.write_text(linkage_to_newick(layout["col_link"],
                                           layout["col_labels"]), encoding="utf-8")
            tree_files["feature_tree"] = f.name
        if not tree_files:
            print("      (Newick needs --cluster and >=2 objects; skipped)",
                  file=sys.stderr)

    # --- manifest ------------------------------------------------------------
    manifest = {
        "mpph_version": __version__,
        "command": "mpph " + " ".join(sys.argv[1:]),
        "python": sys.version.split()[0],
        "platform": sys.platform,
        "source": ("user:" + str(args.user)) if from_user else label_src,
        "mode": mode,
        "metric": metric,
        "match": None if (from_user or args.codes) else match,
        "clustered": bool(args.cluster),
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "kegg_release": (kegg_release(session, cache_dir) if not from_user
                         else "n/a (user data)"),
        "n_organisms_matched": len(organisms),
        "n_organisms": int(df.shape[0]),
        "n_features_before_filtering": n_features_raw,
        "n_features": int(df.shape[1]),
        "organisms": list(df.index),
        "excluded_organisms": (
            [{"organism": o, "reason": "zero_features_after_filtering"}
             for o in excluded_empty] + failed),
        "filters": {
            "mode_metric": metric,
            "top_category": None if (mode != "presence" or args.all_categories)
            else args.top_category,
            "keep_overview": args.keep_overview,
            "min_prevalence": args.min_prevalence,
            "max_prevalence": args.max_prevalence,
            "drop_core": args.drop_core,
            "present_threshold": args.present_threshold,
            "module_types": "all" if args.all_modules else "Pathway",
        },
        "outputs": {
            "matrix": matrix_csv.name,
            "features": features_csv.name,
            "qc": qc_csv.name,
            "figures": figures,
            **ordered_outputs,
            **tree_files,
        },
    }
    (outdir / f"{slug}_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"\nDone. Outputs in {outdir}/")
    for p in sorted(outdir.glob(f"{slug}_*")):
        print(f"  - {p.name}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="mpph",
        description="Metabolic Pathway Presence Heatmap from KEGG (or your MAGs).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("taxon", nargs="?",
                   help="Taxon name to query, e.g. 'Vibrio'. Omit if using "
                        "--codes or --user.")
    src = p.add_argument_group("input source")
    src.add_argument("--match", default="word",
                     choices=["word", "prefix", "exact", "regex"],
                     help="How the taxon matches organism names.")
    src.add_argument("--exact", action="store_true",
                     help="Deprecated alias for --match prefix.")
    src.add_argument("--codes", metavar="FILE",
                     help="File of KEGG organism codes (one per line) to use "
                          "instead of a taxon search.")
    src.add_argument("--user", metavar="PATH",
                     help="Your own KO annotations (a directory of per-MAG KO "
                          "lists, or a KO table). Implies --completeness.")

    ana = p.add_argument_group("analysis")
    ana.add_argument("--completeness", action="store_true",
                     help="Score KEGG *module completeness* (0..1) instead of "
                          "pathway presence/absence.")
    ana.add_argument("--all-modules", action="store_true",
                     help="In completeness mode, score all module types; by "
                          "default only 'Pathway' modules are used (signature "
                          "and reaction modules are excluded).")
    ana.add_argument("--cluster", action="store_true",
                     help="UPGMA-cluster rows/columns and draw dendrograms.")
    ana.add_argument("--metric", default=None,
                     choices=["euclidean", "jaccard", "dice", "hamming",
                              "braycurtis", "cosine"],
                     help="Clustering distance (default: jaccard for presence, "
                          "euclidean for completeness).")
    ana.add_argument("--newick", action="store_true",
                     help="Export the organism and feature trees as Newick.")

    filt = p.add_argument_group("feature filters")
    filt.add_argument("--min-prevalence", type=float, default=0.0,
                      help="Keep features present in >= this fraction of organisms.")
    filt.add_argument("--max-prevalence", type=float, default=1.0,
                      help="Keep features present in <= this fraction of organisms.")
    filt.add_argument("--drop-core", action="store_true",
                      help="Drop features present in ALL organisms (uninformative).")
    filt.add_argument("--keep-overview", action="store_true",
                      help="Keep KEGG 'Global and overview maps' (presence mode).")
    filt.add_argument("--top-category", default="Metabolism",
                      help="In presence mode, keep only pathways under this "
                           "BRITE top-level category.")
    filt.add_argument("--all-categories", action="store_true",
                      help="Keep pathways from all top-level categories, not "
                           "just --top-category.")
    filt.add_argument("--present-threshold", type=float, default=1e-9,
                      help="A feature counts as 'present' for prevalence when "
                           "its value exceeds this (module completeness > 0).")
    filt.add_argument("--keep-empty", action="store_true",
                      help="Keep organisms with zero annotated features (by "
                           "default such genomes are dropped, as an all-absent "
                           "row is unannotated, not truly featureless).")

    out = p.add_argument_group("output")
    out.add_argument("--outdir", default="output", help="Directory for outputs.")
    out.add_argument("--format", nargs="+", default=["pdf"],
                     choices=["pdf", "png", "svg"], help="Figure format(s).")

    cache = p.add_argument_group("cache")
    cache.add_argument("--cache-dir", default=DEFAULT_CACHE_STR,
                       help="Directory for cached KEGG responses.")
    cache.add_argument("--no-cache", action="store_true",
                       help="Disable on-disk caching.")
    cache.add_argument("--refresh", action="store_true",
                       help="Ignore the cache and re-fetch from KEGG.")

    p.add_argument("--version", action="version",
                   version=f"%(prog)s {__version__}")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    sources = [bool(args.taxon), bool(args.codes), bool(args.user)]
    if sum(sources) == 0:
        print("Provide exactly one of: a taxon, --codes FILE, or --user PATH.",
              file=sys.stderr)
        return 2
    if sum(sources) > 1:
        print("Choose only one input source (taxon / --codes / --user).",
              file=sys.stderr)
        return 2
    try:
        return run(args)
    except requests.RequestException as exc:
        print(f"KEGG request failed: {exc}", file=sys.stderr)
        return 1
    except (ValueError, FileNotFoundError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
