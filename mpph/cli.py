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


def _resolve_organisms(args, session, cache_dir):
    """Return ``[(code_or_sample, name), ...]`` from taxon / codes."""
    genomes = list_genomes(session, cache_dir, refresh=args.refresh)
    if args.codes:
        organisms = select_by_codes(genomes, read_code_file(args.codes))
    else:
        organisms = select_by_taxon(genomes, args.taxon, args.exact)
    return organisms


def run(args: argparse.Namespace) -> int:
    cache_dir = None if args.no_cache else Path(args.cache_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    session = make_session()

    from_user = args.user is not None
    mode = "completeness" if (args.completeness or from_user) else "presence"
    label_src = args.user if from_user else (args.taxon or args.codes or "taxon")
    slug = re.sub(r"[^0-9A-Za-z]+", "_", Path(str(label_src)).stem).strip("_") or "mpph"

    # --- gather each organism's feature content ------------------------------
    if from_user:
        from .userdata import load_user_kos
        print(f"[1/4] Loading user KO annotations from {args.user} ...", flush=True)
        org_kos = load_user_kos(args.user)
        organisms = [(name, name) for name in org_kos]
        print(f"      loaded {len(org_kos)} sample(s)")
    else:
        print("[1/4] Fetching KEGG genome list ...", flush=True)
        organisms = _resolve_organisms(args, session, cache_dir)
        if not organisms:
            print(f"No KEGG genomes matched {label_src!r}.", file=sys.stderr)
            print("Try a broader name, drop --exact, or check spelling.",
                  file=sys.stderr)
            return 2
        print(f"      matched {len(organisms)} genome(s)")

    # --- build the matrix ----------------------------------------------------
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
                kos[f"{name} ({code})"] = organism_kos(session, code, cache_dir,
                                                       refresh=args.refresh)
                print(f"      KOs {i}/{len(organisms)}", end="\r", flush=True)
            print()
        df, _ = build_completeness_matrix(kos, module_defs)
    else:
        print(f"[2/4] Fetching pathways for {len(organisms)} genome(s) ...",
              flush=True)
        org_pathways = {}
        for i, (code, name) in enumerate(organisms, 1):
            org_pathways[f"{name} ({code})"] = get_pathways(
                session, code, cache_dir, refresh=args.refresh)
            print(f"      [{i}/{len(organisms)}] {name} ({code})", flush=True)
        categories = fetch_pathway_categories(session, cache_dir,
                                              refresh=args.refresh)
        df, feature_names = build_presence_matrix(org_pathways)

    # --- filter --------------------------------------------------------------
    # Order matters: drop aggregate overview maps, THEN organisms left with no
    # real features (e.g. an unannotated genome carrying only overview maps),
    # and only then apply prevalence/core filtering -- otherwise one all-absent
    # genome caps every prevalence below 1.0 and neuters --drop-core.
    print("[3/4] Filtering matrix ...", flush=True)
    if mode == "presence" and not args.keep_overview:
        overview = {m for m, c in categories.items() if c == OVERVIEW_CATEGORY}
        df = df.drop(columns=[c for c in df.columns if c in overview])
    if not args.keep_empty:
        empty = list(df.index[df.sum(axis=1) == 0])
        if empty:
            print(f"      dropping {len(empty)} unannotated organism(s) with no "
                  f"features: {', '.join(empty)}", flush=True)
            df = df.drop(index=empty)
        if df.shape[0] < 1:
            print("All organisms were empty (no KEGG annotation). Nothing to plot.",
                  file=sys.stderr)
            return 3
    df = filter_matrix(df, args.min_prevalence, args.max_prevalence,
                       args.drop_core)
    if df.shape[1] == 0:
        print("No features left after filtering. Loosen the filters.",
              file=sys.stderr)
        return 3
    print(f"      matrix: {df.shape[0]} organisms x {df.shape[1]} features")

    # --- write tables --------------------------------------------------------
    matrix_csv = outdir / f"{slug}_matrix.csv"
    features_csv = outdir / f"{slug}_features.csv"
    df.to_csv(matrix_csv, index_label="organism")
    kept = list(df.columns)
    pd.DataFrame({
        "feature_id": kept,
        "name": [feature_names.get(k, k) if feature_names else k for k in kept],
        "category": [categories.get(k, "Other") for k in kept],
    }).to_csv(features_csv, index=False)

    # --- figure(s) -----------------------------------------------------------
    print("[4/4] Rendering heatmap ...", flush=True)
    kind = "module completeness" if mode == "completeness" else "pathway presence"
    title = f"Metabolic {kind} · {label_src}"
    subtitle = (f"{df.shape[0]} organisms × {df.shape[1]} "
                + ("KEGG modules" if mode == "completeness" else "KEGG pathways")
                + (f"  ·  UPGMA ({args.metric})" if args.cluster else ""))
    figures, layout = [], {}
    for fmt in args.format:
        fig_path = outdir / f"{slug}_heatmap.{fmt}"
        layout = plot_matrix(df, categories, fig_path, title, subtitle,
                             cluster=args.cluster, mode=mode, metric=args.metric)
        figures.append(fig_path.name)

    # --- Newick tree ---------------------------------------------------------
    newick_file = None
    if args.newick:
        if layout.get("row_link") is not None:
            nwk = linkage_to_newick(layout["row_link"], layout["row_labels"])
            newick_file = outdir / f"{slug}_tree.nwk"
            newick_file.write_text(nwk, encoding="utf-8")
        else:
            print("      (Newick needs --cluster and >=3 organisms; skipped)",
                  file=sys.stderr)

    # --- manifest ------------------------------------------------------------
    manifest = {
        "mpph_version": __version__,
        "source": ("user:" + str(args.user)) if from_user else label_src,
        "mode": mode,
        "metric": args.metric,
        "clustered": bool(args.cluster),
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "kegg_release": (kegg_release(session, cache_dir) if not from_user
                         else "n/a (user data)"),
        "n_organisms": int(df.shape[0]),
        "n_features": int(df.shape[1]),
        "organisms": list(df.index),
        "filters": {
            "min_prevalence": args.min_prevalence,
            "max_prevalence": args.max_prevalence,
            "drop_core": args.drop_core,
            "keep_overview": args.keep_overview,
        },
        "outputs": {
            "matrix": matrix_csv.name,
            "features": features_csv.name,
            "figures": figures,
            "newick": newick_file.name if newick_file else None,
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
    src.add_argument("--exact", action="store_true",
                     help="Match names that START WITH the taxon (a full "
                          "'Genus species'); default is whole-word match.")
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
    ana.add_argument("--metric", default="euclidean",
                     choices=["euclidean", "jaccard", "dice", "hamming"],
                     help="Distance metric for clustering.")
    ana.add_argument("--newick", action="store_true",
                     help="Also export the organism tree as Newick (.nwk).")

    filt = p.add_argument_group("feature filters")
    filt.add_argument("--min-prevalence", type=float, default=0.0,
                      help="Keep features present in >= this fraction of organisms.")
    filt.add_argument("--max-prevalence", type=float, default=1.0,
                      help="Keep features present in <= this fraction of organisms.")
    filt.add_argument("--drop-core", action="store_true",
                      help="Drop features present in ALL organisms (uninformative).")
    filt.add_argument("--keep-overview", action="store_true",
                      help="Keep KEGG 'Global and overview maps' (presence mode).")
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
    if not args.taxon and not args.codes and not args.user:
        print("Provide a taxon, --codes FILE, or --user PATH.", file=sys.stderr)
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
