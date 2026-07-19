"""Command-line interface for MPPH."""
from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
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
from .plot import (
    plot_accumulation,
    plot_enrichment,
    plot_gsea_running,
    plot_gsea_summary,
    plot_matrix,
    plot_ordination,
    plot_prevalence,
    plot_volcano,
)
from .tree import linkage_to_newick

DEFAULT_CACHE_STR = str(DEFAULT_CACHE)
SUBCOMMANDS = ("run", "traits", "explain", "compare", "pan", "ordination",
               "community", "report", "validate", "enrich", "gsea")


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
                       args.drop_core, present_threshold=args.present_threshold,
                       prevalence_state=args.prevalence_state,
                       complete_threshold=args.complete_threshold)
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
    strip_label = ("KEGG module functional category" if mode == "completeness"
                   else "KEGG functional category")
    figures, layout = [], {}
    for fmt in args.format:
        fig_path = outdir / f"{slug}_heatmap.{fmt}"
        layout = plot_matrix(df, categories, fig_path, title, subtitle,
                             cluster=args.cluster, mode=mode, metric=metric,
                             strip_label=strip_label)
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
                                           layout["row_link_labels"]), encoding="utf-8")
            tree_files["organism_tree"] = f.name
        if layout.get("col_link") is not None:
            f = outdir / f"{slug}_feature_tree.nwk"
            f.write_text(linkage_to_newick(layout["col_link"],
                                           layout["col_link_labels"]), encoding="utf-8")
            tree_files["feature_tree"] = f.name
        if not tree_files:
            print("      (Newick needs --cluster and >=2 objects; skipped)",
                  file=sys.stderr)

    # --- manifest ------------------------------------------------------------
    manifest = {
        "mpph_version": __version__,
        # the argv actually passed in (sys.argv is wrong for main(argv=...))
        "command": "mpph " + " ".join(getattr(args, "_raw_argv", sys.argv[1:])),
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
            "prevalence_state": args.prevalence_state,
            "complete_threshold": args.complete_threshold,
            "module_types": "all" if args.all_modules else "Pathway",
        },
        "score_semantics": ("module_step_coverage" if mode == "completeness"
                            else "pathway_map_association"),
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
    problems: list[str] = []
    if has_files:
        # Relative annotation_file paths are resolved against the sheet's own
        # directory -- the same convention load_kos_from_sheet's caller uses.
        base = Path(args.samples).parent
        raw = sheet["annotation_file"].astype(str).str.strip()
        blank = raw == ""
        if blank.any():
            problems.append(
                f"{int(blank.sum())} row(s) with an empty annotation_file: "
                f"{sheet.loc[blank, 'sample_id'].tolist()}")
        for sid, af in zip(sheet["sample_id"], raw):
            if af and not (base / af).exists():
                problems.append(f"annotation_file not found for "
                                f"sample_id={sid!r}: {base / af}")
    if problems:
        for p in problems:
            print(f"PROBLEM: {p}", file=sys.stderr)
        return 1
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
    acc = accumulation_curve(matrix, permutations=args.permutations,
                             seed=args.seed)
    acc.to_csv(results / f"{slug}_accumulation.csv", index=False)
    plot_prevalence(classes, results / f"{slug}_prevalence.{args.format}",
                    f"Pan-functional prevalence · {slug}")
    plot_accumulation(acc, results / f"{slug}_accumulation.{args.format}",
                      f"Functional accumulation · {slug}")
    counts = classes["pan_class"].value_counts().to_dict()
    print(f"Pan-functional classes: {counts}")
    print(f"Wrote {slug}_pan_classes.csv, _accumulation.csv, "
          f"_prevalence.{args.format}, _accumulation.{args.format}")
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
    stem = f"{slug}_differential_{args.group_a}_vs_{args.group_b}"
    out.to_csv(results / f"{stem}.csv", index=False)
    plot_volcano(out, results / f"{stem}.{args.format}",
                 f"{args.group_a} vs {args.group_b} · {slug}")
    n_sig = int((out["q_value"] < 0.05).sum())
    print(f"{len(out)} features tested, {n_sig} with q<0.05. "
          f"Wrote {stem}.csv and {stem}.{args.format}")
    return 0


def cmd_ordination(args) -> int:
    from .analysis import pcoa, permanova
    results = Path(args.results)
    slug = _find_slug(results, args.slug)
    matrix = _read_matrix(results, slug)
    coords, explained, diagnostics = pcoa(matrix, metric=args.metric)
    coords.to_csv(results / f"{slug}_pcoa.csv", index_label="organism")
    (results / f"{slug}_pcoa_diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2), encoding="utf-8")
    groups = None
    if args.metadata and args.color:
        groups = _read_metadata(Path(args.metadata))[args.color]
    plot_ordination(coords, explained, results / f"{slug}_pcoa.{args.format}",
                    f"PCoA ({args.metric}) · {slug}", groups)
    pc2 = explained[1] if len(explained) > 1 else 0.0
    print(f"PCoA axis-1/2 variance: {explained[0]:.1%} / {pc2:.1%} "
          f"(negative eigenvalue fraction {diagnostics['negative_fraction']:.1%})")
    if groups is not None:
        res = permanova(matrix, groups, metric=args.metric,
                        permutations=args.permutations, seed=args.seed)
        pd.DataFrame([res]).to_csv(results / f"{slug}_permanova.csv", index=False)
        print(f"PERMANOVA: pseudo-F={res['pseudo_F']:.3f}, p={res['p_value']:.4f}"
              + (f" ({res['n_excluded_missing_label']} sample(s) excluded for a "
                 "missing group label)" if res["n_excluded_missing_label"] else ""))
    print(f"Wrote {slug}_pcoa.csv, _pcoa_diagnostics.json and _pcoa.{args.format}")
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
        mark = ("OK " if step.satisfied is True else
                "?? " if step.satisfied is None else "-- ")
        print(f"  {mark}step {i}: {step.expression}")
        if step.matched_kos:
            print(f"       matched: {', '.join(step.matched_kos)}")
        if step.satisfied is False and step.missing_kos:
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
                    metric="euclidean", value_label="trait completeness",
                    strip_label="trait category")
    print(f"Scored {matrix.shape[1]} traits for {matrix.shape[0]} organisms. "
          f"Wrote to {outdir}/")
    return 0


def cmd_enrich(args) -> int:
    from .enrichment import (
        fetch_ko_module_membership,
        fetch_ko_pathway_membership,
        fetch_pathway_names,
        hypergeometric_enrichment,
        invert_membership,
        load_gene_go_map,
        load_go_names,
        read_id_list,
    )

    if args.ontology == "go" and not args.gene_go_map:
        print("Error: --ontology go requires --gene-go-map FILE (KEGG has no "
              "GO annotations; supply your own gene-to-GO mapping).",
              file=sys.stderr)
        return 2

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    study = read_id_list(args.study)

    if args.ontology == "go":
        if args.background_organism:
            print("Warning: --background-organism has no meaning for "
                  "--ontology go (GO mapping is not organism-specific); "
                  "use --background instead.", file=sys.stderr)
        item_to_cats = load_gene_go_map(args.gene_go_map, args.go_map_format)
        category_names = load_go_names(args.go_names) if args.go_names else {}
        if not args.background:
            print("Error: --ontology go requires --background FILE.",
                  file=sys.stderr)
            return 2
        background = read_id_list(args.background)
    else:
        cache_dir = None if args.no_cache else Path(args.cache_dir)
        session = make_session()
        if args.ontology == "kegg-pathway":
            item_to_cats = fetch_ko_pathway_membership(session, cache_dir,
                                                       refresh=args.refresh)
            category_names = fetch_pathway_names(session, cache_dir,
                                                 refresh=args.refresh)
        else:  # kegg-module
            item_to_cats = fetch_ko_module_membership(session, cache_dir,
                                                       refresh=args.refresh)
            category_names = list_modules(session, cache_dir, refresh=args.refresh)
        if args.background_organism:
            background = organism_kos(session, args.background_organism,
                                      cache_dir, refresh=args.refresh)
        elif args.background:
            background = read_id_list(args.background)
        else:
            print("Error: provide --background FILE or --background-organism CODE.",
                  file=sys.stderr)
            return 2

    category_to_items = invert_membership(item_to_cats)
    results, stats = hypergeometric_enrichment(
        study, background, category_to_items, category_names,
        min_category_size=args.min_category_size)

    dest_csv = outdir / f"{args.label}_enrichment.csv"
    results.to_csv(dest_csv, index=False)
    print(f"Study: {stats['n_study_used']}/{stats['n_study_input']} gene(s) "
          f"annotated & in background ({stats['n_study_not_in_background']} not "
          f"in background). Background: {stats['n_background_used']}/"
          f"{stats['n_background_input']} annotated. Tested "
          f"{stats['n_categories_tested']} categories (>= "
          f"{args.min_category_size} background members).")
    n_sig = int((results["q_value"] < args.alpha).sum()) if len(results) else 0
    print(f"{n_sig} categories significant at q<{args.alpha}. Wrote {dest_csv.name}")

    if len(results):
        dest_fig = outdir / f"{args.label}_enrichment.{args.format}"
        plot_enrichment(results, dest_fig,
                        f"Enrichment ({args.ontology}) · {args.label}",
                        top_n=args.top_n, alpha=args.alpha)
        print(f"Wrote {dest_fig.name}")
    else:
        print("No category met --min-category-size; skipping the plot.",
              file=sys.stderr)
    return 0


def cmd_gsea(args) -> int:
    from .enrichment import (
        fetch_ko_module_membership,
        fetch_ko_pathway_membership,
        fetch_pathway_names,
        invert_membership,
        load_gene_go_map,
        load_go_names,
    )
    from .gsea import (
        gsea_analysis,
        load_expression_matrix,
        load_ranked_list,
        rank_from_expression,
    )

    if args.ontology == "go" and not args.gene_go_map:
        print("Error: --ontology go requires --gene-go-map FILE (KEGG has no "
              "GO annotations; supply your own gene-to-GO mapping).",
              file=sys.stderr)
        return 2

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.ranked_list:
        ranked = load_ranked_list(args.ranked_list)
    else:
        expr = load_expression_matrix(args.expression)
        groups = _read_metadata(Path(args.metadata))[args.group_column]
        ranked = rank_from_expression(expr, groups, args.group_a, args.group_b,
                                      metric=args.rank_metric)
        ranked.to_csv(outdir / f"{args.label}_ranked_list.tsv", sep="\t",
                     header=False)

    if args.ontology == "go":
        item_to_cats = load_gene_go_map(args.gene_go_map, args.go_map_format)
        category_names = load_go_names(args.go_names) if args.go_names else {}
    else:
        cache_dir = None if args.no_cache else Path(args.cache_dir)
        session = make_session()
        if args.ontology == "kegg-pathway":
            item_to_cats = fetch_ko_pathway_membership(session, cache_dir,
                                                       refresh=args.refresh)
            category_names = fetch_pathway_names(session, cache_dir,
                                                 refresh=args.refresh)
        else:  # kegg-module
            item_to_cats = fetch_ko_module_membership(session, cache_dir,
                                                       refresh=args.refresh)
            category_names = list_modules(session, cache_dir, refresh=args.refresh)

    category_to_items = invert_membership(item_to_cats)
    results, running_sums, ranked_genes = gsea_analysis(
        ranked, category_to_items, category_names,
        weight=args.weight, min_size=args.min_size, max_size=args.max_size,
        permutations=args.permutations, seed=args.seed)

    dest_csv = outdir / f"{args.label}_gsea.csv"
    results.to_csv(dest_csv, index=False)
    n_sig = int((results["q_value"] < args.alpha).sum()) if len(results) else 0
    print(f"Ranked {len(ranked)} genes. Tested {len(results)} categories "
          f"(size {args.min_size}-{args.max_size} genes). "
          f"{n_sig} significant at q<{args.alpha}. Wrote {dest_csv.name}")

    if len(results):
        dest_summary = outdir / f"{args.label}_gsea_summary.{args.format}"
        plot_gsea_summary(results, dest_summary,
                          f"GSEA ({args.ontology}) · {args.label}",
                          top_n=args.top_n, alpha=args.alpha)
        print(f"Wrote {dest_summary.name}")

        top = results.iloc[0]
        top_members = category_to_items.get(top["category_id"], set())
        hits = np.fromiter((g in top_members for g in ranked_genes), dtype=bool,
                          count=len(ranked_genes))
        dest_run = outdir / f"{args.label}_gsea_top.{args.format}"
        plot_gsea_running(
            ranked.to_numpy(dtype=float), running_sums[top["category_id"]], hits,
            dest_run, f"{top['category_id']}  {top['category_name']}",
            f"NES={top['NES']:.2f}  q={top['q_value']:.3g}  "
            f"({top['leading_edge_size']}/{top['size']} leading-edge genes)")
        print(f"Wrote {dest_run.name}")
    else:
        print("No category met the size filters; skipping plots.",
              file=sys.stderr)
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
                     choices=["auto", "ko-list", "long", "eggnog"],
                     help="Format of --user input. 'long' forces sample<TAB>KO "
                          "table parsing for an ambiguous two-column file "
                          "where 'auto' would otherwise treat each row as its "
                          "own sample.")
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
    filt.add_argument("--prevalence-state", default="any",
                      choices=["any", "complete"],
                      help="Count a feature toward prevalence when it is "
                           "detectable ('any') or fully complete ('complete').")
    filt.add_argument("--complete-threshold", type=float, default=1.0,
                      help="Completeness at/above which a module counts as "
                           "'complete' for --prevalence-state complete.")
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
                   choices=["auto", "ko-list", "long", "eggnog"])
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
    cmp.add_argument("--format", default="png", choices=["pdf", "png", "svg"])

    pan = sub.add_parser("pan", help="Core/soft-core/shell/cloud classification.")
    pan.add_argument("--results", required=True)
    pan.add_argument("--slug")
    pan.add_argument("--core", type=float, default=1.0)
    pan.add_argument("--soft-core", type=float, default=0.95)
    pan.add_argument("--shell", type=float, default=0.15)
    pan.add_argument("--permutations", type=int, default=100)
    pan.add_argument("--seed", type=int, default=0)
    pan.add_argument("--format", default="png", choices=["pdf", "png", "svg"])

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

    enr = sub.add_parser(
        "enrich", help="KEGG pathway/module or GO over-representation (ORA).")
    enr.add_argument("--study", required=True, metavar="FILE",
                     help="Study-set gene/KO ids, one per line.")
    enr.add_argument("--background", metavar="FILE",
                     help="Background/universe gene or KO ids, one per line "
                          "(required for --ontology go).")
    enr.add_argument("--background-organism", metavar="CODE",
                     help="KEGG organism code; use its full KO complement as "
                          "the background (kegg-pathway/kegg-module only).")
    enr.add_argument("--ontology", required=True,
                     choices=["kegg-pathway", "kegg-module", "go"],
                     help="Category source to test the study set against.")
    enr.add_argument("--gene-go-map", metavar="FILE",
                     help="Required for --ontology go: a gene-to-GO mapping "
                          "(long table, or an eggNOG-mapper .annotations file "
                          "-- KEGG itself has no GO annotations).")
    enr.add_argument("--go-map-format", default="auto",
                     choices=["auto", "eggnog"])
    enr.add_argument("--go-names", metavar="FILE",
                     help="Optional GO id -> name table for readable labels.")
    enr.add_argument("--min-category-size", type=int, default=2,
                     help="Skip categories with fewer background members.")
    enr.add_argument("--top-n", type=int, default=20,
                     help="Categories shown in the plot.")
    enr.add_argument("--alpha", type=float, default=0.05,
                     help="Significance threshold drawn on the plot.")
    enr.add_argument("--label", default="enrichment",
                     help="Prefix for output file names.")
    enr.add_argument("--outdir", default="output")
    enr.add_argument("--format", default="png", choices=["pdf", "png", "svg"])
    _add_cache_args(enr)

    gse = sub.add_parser(
        "gsea", help="Rank-based (GSEA-style) enrichment from expression or a "
                     "pre-ranked gene list.")
    rank_src = gse.add_mutually_exclusive_group(required=True)
    rank_src.add_argument("--ranked-list", metavar="FILE",
                          help="Pre-ranked gene<TAB>score file (e.g. your own "
                               "DESeq2/edgeR/limma statistic).")
    rank_src.add_argument("--expression", metavar="FILE",
                          help="gene<TAB>sample... expression matrix; ranks "
                               "genes automatically with --metadata/--group-*. "
                               "Must already be normalized for library size "
                               "(CPM/TPM/FPKM or DESeq2/edgeR size factors) -- "
                               "raw read counts are not comparable across "
                               "samples sequenced to different depths.")
    gse.add_argument("--metadata", metavar="FILE",
                     help="Sample metadata (needs sample_id + --group-column); "
                          "required with --expression.")
    gse.add_argument("--group-column", help="Metadata column with group labels.")
    gse.add_argument("--group-a", help="First group (ranked toward the top).")
    gse.add_argument("--group-b", help="Second group (ranked toward the bottom).")
    gse.add_argument("--rank-metric", default="signal2noise",
                     choices=["signal2noise", "log2fc"],
                     help="Ranking statistic when using --expression.")
    gse.add_argument("--ontology", required=True,
                     choices=["kegg-pathway", "kegg-module", "go"],
                     help="Category source to test the ranking against.")
    gse.add_argument("--gene-go-map", metavar="FILE",
                     help="Required for --ontology go (see `enrich --help`).")
    gse.add_argument("--go-map-format", default="auto",
                     choices=["auto", "eggnog"])
    gse.add_argument("--go-names", metavar="FILE",
                     help="Optional GO id -> name table for readable labels.")
    gse.add_argument("--min-size", type=int, default=15,
                     help="Skip categories with fewer members in the ranking.")
    gse.add_argument("--max-size", type=int, default=500,
                     help="Skip categories with more members in the ranking.")
    gse.add_argument("--weight", type=float, default=1.0,
                     help="Score exponent in the running-sum statistic (0 = "
                          "unweighted KS; 1 = standard GSEA weighting).")
    gse.add_argument("--permutations", type=int, default=1000,
                     help="Gene-set permutations per distinct category size.")
    gse.add_argument("--seed", type=int, default=0)
    gse.add_argument("--top-n", type=int, default=20,
                     help="Categories shown in the summary plot.")
    gse.add_argument("--alpha", type=float, default=0.05)
    gse.add_argument("--label", default="gsea", help="Prefix for output files.")
    gse.add_argument("--outdir", default="output")
    gse.add_argument("--format", default="png", choices=["pdf", "png", "svg"])
    _add_cache_args(gse)
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
    args._raw_argv = raw  # faithful provenance, even when called as main(argv=...)

    handlers = {
        "run": run, "traits": cmd_traits, "explain": cmd_explain,
        "compare": cmd_compare, "pan": cmd_pan, "ordination": cmd_ordination,
        "community": cmd_community, "report": cmd_report, "validate": cmd_validate,
        "enrich": cmd_enrich, "gsea": cmd_gsea,
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
    if args.command == "gsea" and args.expression and not (
            args.metadata and args.group_column and args.group_a and args.group_b):
        print("Error: --expression needs --metadata, --group-column, "
              "--group-a and --group-b to rank genes.", file=sys.stderr)
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
