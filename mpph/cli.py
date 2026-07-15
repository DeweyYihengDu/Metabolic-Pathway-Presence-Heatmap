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
from .modules import evaluate_module
from .pathways import OVERVIEW_CATEGORY, fetch_pathway_categories, get_pathways
from .plot import plot_matrix, plot_ordination
from .tree import linkage_to_newick

DEFAULT_CACHE_STR = str(DEFAULT_CACHE)
SUBCOMMANDS = ("run", "traits", "explain", "compare", "pan", "ordination",
               "community", "report", "validate")


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
        org_kos = load_user_kos(args.user, args.input_format)
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

    if args.dry_run:
        print("Matched organisms (dry run, no data fetched):")
        for code, name in organisms:
            print(f"  {code}\t{name}")
        return 0

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


# --------------------------------------------------------------------------- #
# Helpers shared by KEGG-touching subcommands
# --------------------------------------------------------------------------- #
def _load_source_kos(args, session, cache_dir) -> dict[str, set[str]]:
    """Return ``{label: KO set}`` from a taxon / --codes / --user source."""
    if args.user:
        from .userdata import load_user_kos
        return load_user_kos(args.user, args.input_format)
    genomes = list_genomes(session, cache_dir, refresh=args.refresh)
    if args.codes:
        organisms = select_by_codes(genomes, read_code_file(args.codes))
    else:
        match = "prefix" if args.exact else args.match
        organisms = select_by_taxon(genomes, args.taxon, match)
    if not organisms:
        raise ValueError(f"no organisms matched {args.taxon or args.codes!r}")
    return {f"{name} ({code})": organism_kos(session, code, cache_dir,
                                             refresh=args.refresh)
            for code, name in organisms}


def _load_pathway_module_defs(args, session, cache_dir):
    modules = list_modules(session, cache_dir, refresh=args.refresh)
    defs = fetch_module_definitions(session, sorted(modules), cache_dir,
                                    refresh=args.refresh)
    if not getattr(args, "all_modules", False):
        defs = {m: v for m, v in defs.items() if v[2] == "Pathway"}
    return defs


def _read_matrix(results: Path, slug: str) -> pd.DataFrame:
    df = pd.read_csv(results / f"{slug}_matrix.csv", index_col=0)
    df.columns = df.columns.astype(str)
    return df


def _find_slug(results: Path, slug: str | None) -> str:
    if slug:
        return slug
    hits = sorted(results.glob("*_matrix.csv"))
    if not hits:
        raise ValueError(f"no *_matrix.csv found in {results}; pass --slug")
    return hits[0].name[: -len("_matrix.csv")]


def _read_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep=None, engine="python", dtype=str)
    key = "sample_id" if "sample_id" in meta.columns else meta.columns[0]
    return meta.set_index(key)


# --------------------------------------------------------------------------- #
# Subcommands (post-processing an output folder or recomputing from a source)
# --------------------------------------------------------------------------- #
def cmd_report(args) -> int:
    from .report import build_report
    results = Path(args.results)
    slug = _find_slug(results, args.slug)
    out = build_report(results, slug)
    print(f"Wrote {out}")
    return 0


def cmd_validate(args) -> int:
    from .samplesheet import read_sample_sheet
    sheet = read_sample_sheet(args.samples)
    n = len(sheet)
    has_files = "annotation_file" in sheet.columns
    print(f"OK: {n} samples, unique ids"
          + (", annotation_file present" if has_files else
             " (no annotation_file column)"))
    return 0


def cmd_pan(args) -> int:
    from .analysis import accumulation_curve, pan_classify
    results = Path(args.results)
    slug = _find_slug(results, args.slug)
    matrix = _read_matrix(results, slug)
    classes = pan_classify(matrix, core=args.core, soft_core=args.soft_core,
                           shell=args.shell)
    classes.to_csv(results / f"{slug}_pan_classes.csv", index=False)
    accumulation_curve(matrix, permutations=args.permutations, seed=args.seed) \
        .to_csv(results / f"{slug}_accumulation.csv", index=False)
    counts = classes["pan_class"].value_counts().to_dict()
    print(f"Pan-functional classes: {counts}")
    print(f"Wrote {slug}_pan_classes.csv and {slug}_accumulation.csv")
    return 0


def cmd_compare(args) -> int:
    from .analysis import differential_features
    results = Path(args.results)
    slug = _find_slug(results, args.slug)
    matrix = _read_matrix(results, slug)
    meta = _read_metadata(Path(args.metadata))
    groups = meta[args.group_column]
    manifest_path = results / f"{slug}_manifest.json"
    continuous = (json.loads(manifest_path.read_text("utf-8")).get("mode")
                  == "completeness") if manifest_path.exists() else False
    out = differential_features(matrix, groups, args.group_a, args.group_b,
                                continuous=continuous)
    dest = results / f"{slug}_differential_{args.group_a}_vs_{args.group_b}.csv"
    out.to_csv(dest, index=False)
    n_sig = int((out["q_value"] < 0.05).sum())
    print(f"{len(out)} features tested, {n_sig} with q<0.05. Wrote {dest.name}")
    return 0


def cmd_ordination(args) -> int:
    from .analysis import pcoa, permanova
    results = Path(args.results)
    slug = _find_slug(results, args.slug)
    matrix = _read_matrix(results, slug)
    coords, explained = pcoa(matrix, metric=args.metric)
    coords.to_csv(results / f"{slug}_pcoa.csv", index_label="organism")
    groups = None
    if args.metadata and args.color:
        groups = _read_metadata(Path(args.metadata))[args.color]
    plot_ordination(coords, explained, results / f"{slug}_pcoa.{args.format}",
                    f"PCoA ({args.metric}) · {slug}", groups)
    print(f"PCoA axis-1/2 variance: {explained[0]:.1%} / "
          f"{explained[1] if len(explained) > 1 else 0:.1%}")
    if groups is not None:
        res = permanova(matrix, groups, metric=args.metric,
                        permutations=args.permutations, seed=args.seed)
        print(f"PERMANOVA: pseudo-F={res['pseudo_F']:.3f}, p={res['p_value']:.4f}")
    print(f"Wrote {slug}_pcoa.csv and {slug}_pcoa.{args.format}")
    return 0


def cmd_explain(args) -> int:
    cache_dir = None if args.no_cache else Path(args.cache_dir)
    session = make_session()
    defs = _load_pathway_module_defs(args, session, cache_dir)
    if args.module not in defs:
        # fetch this single module even if not a Pathway module
        extra = fetch_module_definitions(session, [args.module], cache_dir,
                                         refresh=args.refresh)
        defs = {**defs, **extra}
    if args.module not in defs:
        print(f"Module {args.module} not found.", file=sys.stderr)
        return 2
    kos = _load_source_kos(args, session, cache_dir)
    label = args.organism or next(iter(kos))
    match = [k for k in kos if label in k or k == label]
    if not match:
        print(f"Organism {label!r} not among: {list(kos)}", file=sys.stderr)
        return 2
    ko_set = kos[match[0]]
    ev = evaluate_module(args.module, defs[args.module][0], ko_set, defs)
    print(f"Sample: {match[0]}")
    print(f"Module: {ev.module_id}  ({defs[args.module][1]})")
    print(f"State: {ev.state}   coverage: {ev.n_satisfied}/{ev.n_steps} "
          f"= {ev.score:.3f}   parser: {ev.parser_status}")
    for i, step in enumerate(ev.steps, 1):
        mark = "OK " if step.satisfied else "-- "
        print(f"  {mark}step {i}: {step.expression}")
        if step.matched_kos:
            print(f"       matched: {', '.join(step.matched_kos)}")
        if not step.satisfied and step.missing_kos:
            print(f"       missing: {', '.join(step.missing_kos)}")
    if ev.unresolved_references:
        print(f"  unresolved module refs: {', '.join(ev.unresolved_references)}")
    return 0


def cmd_community(args) -> int:
    from .analysis import pairwise_complementarity
    cache_dir = None if args.no_cache else Path(args.cache_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    session = make_session()
    org_kos = _load_source_kos(args, session, cache_dir)
    defs = _load_pathway_module_defs(args, session, cache_dir)
    out = pairwise_complementarity(org_kos, defs,
                                   max_combination_size=args.max_combination_size,
                                   only_newly_completed=not args.all_pairs)
    dest = outdir / "community_complementarity.csv"
    out.to_csv(dest, index=False)
    print(f"{len(out)} complementarity events. Wrote {dest}")
    return 0


def cmd_traits(args) -> int:
    from .traits import load_trait_panel, score_traits
    cache_dir = None if args.no_cache else Path(args.cache_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    session = make_session()
    org_kos = _load_source_kos(args, session, cache_dir)
    panel = load_trait_panel(args.panel)
    matrix, names, categories = score_traits(org_kos, panel)
    slug = re.sub(r"[^0-9A-Za-z]+", "_", Path(args.panel).stem) or "traits"
    matrix.to_csv(outdir / f"{slug}_matrix.csv", index_label="organism")
    title = f"Metabolic traits · {slug}"
    subtitle = f"{matrix.shape[0]} organisms × {matrix.shape[1]} traits"
    for fmt in args.format:
        plot_matrix(matrix, categories, outdir / f"{slug}_heatmap.{fmt}",
                    title, subtitle, cluster=args.cluster, mode="completeness",
                    metric="euclidean")
    print(f"Scored {matrix.shape[1]} traits for {matrix.shape[0]} organisms. "
          f"Wrote to {outdir}/")
    return 0


# --------------------------------------------------------------------------- #
# Argument parser (subcommands; `mpph <taxon> ...` is a shortcut for `run`)
# --------------------------------------------------------------------------- #
def _add_run_arguments(p: argparse.ArgumentParser) -> None:
    p.add_argument("taxon", nargs="?",
                   help="Taxon name to query, e.g. 'Vibrio'. Omit for --codes/--user.")
    src = p.add_argument_group("input source")
    src.add_argument("--match", default="word",
                     choices=["word", "prefix", "exact", "regex"],
                     help="How the taxon matches organism names.")
    src.add_argument("--exact", action="store_true",
                     help="Deprecated alias for --match prefix.")
    src.add_argument("--codes", metavar="FILE",
                     help="File of KEGG organism codes (one per line).")
    src.add_argument("--user", metavar="PATH",
                     help="Your own KO annotations. Implies --completeness.")
    src.add_argument("--input-format", default="auto",
                     choices=["auto", "ko-list", "eggnog"],
                     help="Format of --user input.")
    ana = p.add_argument_group("analysis")
    ana.add_argument("--completeness", action="store_true",
                     help="Score KEGG module completeness (0..1).")
    ana.add_argument("--all-modules", action="store_true",
                     help="Score all module types (default: Pathway only).")
    ana.add_argument("--cluster", action="store_true",
                     help="UPGMA-cluster rows/columns and draw dendrograms.")
    ana.add_argument("--metric", default=None,
                     choices=["euclidean", "jaccard", "dice", "hamming",
                              "braycurtis", "cosine"],
                     help="Clustering distance (default per mode).")
    ana.add_argument("--newick", action="store_true",
                     help="Export organism and feature trees as Newick.")
    ana.add_argument("--dry-run", action="store_true",
                     help="List matched organisms and exit.")
    filt = p.add_argument_group("feature filters")
    filt.add_argument("--min-prevalence", type=float, default=0.0)
    filt.add_argument("--max-prevalence", type=float, default=1.0)
    filt.add_argument("--drop-core", action="store_true")
    filt.add_argument("--keep-overview", action="store_true")
    filt.add_argument("--top-category", default="Metabolism")
    filt.add_argument("--all-categories", action="store_true")
    filt.add_argument("--present-threshold", type=float, default=1e-9)
    filt.add_argument("--keep-empty", action="store_true")
    out = p.add_argument_group("output")
    out.add_argument("--outdir", default="output")
    out.add_argument("--format", nargs="+", default=["pdf"],
                     choices=["pdf", "png", "svg"])
    _add_cache_args(p)


def _add_source_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("taxon", nargs="?", help="Taxon name to query.")
    p.add_argument("--match", default="word",
                   choices=["word", "prefix", "exact", "regex"])
    p.add_argument("--exact", action="store_true")
    p.add_argument("--codes", metavar="FILE")
    p.add_argument("--user", metavar="PATH")
    p.add_argument("--input-format", default="auto",
                   choices=["auto", "ko-list", "eggnog"])
    p.add_argument("--all-modules", action="store_true")


def _add_cache_args(p: argparse.ArgumentParser) -> None:
    cache = p.add_argument_group("cache")
    cache.add_argument("--cache-dir", default=DEFAULT_CACHE_STR)
    cache.add_argument("--no-cache", action="store_true")
    cache.add_argument("--refresh", action="store_true")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="mpph",
        description="KEGG functional-profile analysis: heatmaps, module "
                    "completeness, traits, comparison, community, reports.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="command")

    _add_run_arguments(sub.add_parser("run", help="Build a heatmap (default)."))

    tr = sub.add_parser("traits", help="Score a metabolic-trait panel.")
    _add_source_args(tr)
    tr.add_argument("--panel", default="biogeochemistry",
                    help="Trait panel name (built-in) or path to JSON/YAML.")
    tr.add_argument("--cluster", action="store_true")
    tr.add_argument("--outdir", default="output")
    tr.add_argument("--format", nargs="+", default=["png"],
                    choices=["pdf", "png", "svg"])
    _add_cache_args(tr)

    ex = sub.add_parser("explain", help="Explain one module's evidence.")
    _add_source_args(ex)
    ex.add_argument("--module", required=True, help="Module id, e.g. M00002.")
    ex.add_argument("--organism", help="Organism label substring (default: first).")
    _add_cache_args(ex)

    cmp = sub.add_parser("compare", help="Differential features between groups.")
    cmp.add_argument("--results", required=True)
    cmp.add_argument("--slug")
    cmp.add_argument("--metadata", required=True)
    cmp.add_argument("--group-column", required=True)
    cmp.add_argument("--group-a", required=True)
    cmp.add_argument("--group-b", required=True)

    pan = sub.add_parser("pan", help="Core/soft-core/shell/cloud classification.")
    pan.add_argument("--results", required=True)
    pan.add_argument("--slug")
    pan.add_argument("--core", type=float, default=1.0)
    pan.add_argument("--soft-core", type=float, default=0.95)
    pan.add_argument("--shell", type=float, default=0.15)
    pan.add_argument("--permutations", type=int, default=100)
    pan.add_argument("--seed", type=int, default=0)

    ordn = sub.add_parser("ordination", help="PCoA + optional PERMANOVA.")
    ordn.add_argument("--results", required=True)
    ordn.add_argument("--slug")
    ordn.add_argument("--metric", default="braycurtis")
    ordn.add_argument("--metadata")
    ordn.add_argument("--color")
    ordn.add_argument("--permutations", type=int, default=999)
    ordn.add_argument("--seed", type=int, default=0)
    ordn.add_argument("--format", default="png", choices=["pdf", "png", "svg"])

    com = sub.add_parser("community", help="Pairwise metabolic complementarity.")
    _add_source_args(com)
    com.add_argument("--max-combination-size", type=int, default=2)
    com.add_argument("--all-pairs", action="store_true",
                     help="Report all pairs, not only newly completed modules.")
    com.add_argument("--outdir", default="output")
    _add_cache_args(com)

    rep = sub.add_parser("report", help="Build a self-contained HTML report.")
    rep.add_argument("--results", required=True)
    rep.add_argument("--slug")

    val = sub.add_parser("validate", help="Validate a sample sheet.")
    val.add_argument("--samples", required=True)
    return p


def _validate_source(args) -> int | None:
    sources = [bool(args.taxon), bool(args.codes), bool(args.user)]
    if sum(sources) == 0:
        print("Provide exactly one of: a taxon, --codes FILE, or --user PATH.",
              file=sys.stderr)
        return 2
    if sum(sources) > 1:
        print("Choose only one input source (taxon / --codes / --user).",
              file=sys.stderr)
        return 2
    return None


def main(argv: list[str] | None = None) -> int:
    raw = list(sys.argv[1:] if argv is None else argv)
    # Backward-compatible shortcut: `mpph <taxon> ...` -> `mpph run <taxon> ...`
    if raw and not raw[0].startswith("-") and raw[0] not in SUBCOMMANDS:
        raw = ["run", *raw]
    args = build_parser().parse_args(raw)
    if not args.command:
        build_parser().print_help(sys.stderr)
        return 2

    handlers = {
        "run": run, "traits": cmd_traits, "explain": cmd_explain,
        "compare": cmd_compare, "pan": cmd_pan, "ordination": cmd_ordination,
        "community": cmd_community, "report": cmd_report, "validate": cmd_validate,
    }
    if args.command in ("run", "traits", "explain", "community"):
        err = _validate_source(args)
        if err:
            return err
    if args.command == "run" and not (
            0.0 <= args.min_prevalence <= args.max_prevalence <= 1.0):
        print("Error: require 0 <= --min-prevalence <= --max-prevalence <= 1.",
              file=sys.stderr)
        return 2
    try:
        return handlers[args.command](args)
    except requests.RequestException as exc:
        print(f"KEGG request failed: {exc}", file=sys.stderr)
        return 1
    except (ValueError, FileNotFoundError, KeyError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
