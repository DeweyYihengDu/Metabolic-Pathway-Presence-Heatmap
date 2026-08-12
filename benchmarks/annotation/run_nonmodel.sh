#!/usr/bin/env bash
# Extends run_all.sh's three model organisms with the harder cases:
#
#   1. Non-model, environmentally-derived genomes KEGG never hand-curated
#      (`vbs` Verrucomicrobia bacterium S94, `lbac` Lentisphaerae bacterium
#      WC36 -- both PVC superphylum, both added to KEGG only recently, both
#      with placeholder "bacterium <strain>" names rather than species
#      names). These still have real KEGG KO assignments, so precision/
#      recall against ground truth is still measured.
#
#   2. Real marine metagenome-assembled genomes (MAGs, Prodigal gene calls)
#      -- the actual target use case for `mpph annotate`: fragmentary,
#      incomplete, no reference annotation whatsoever. No ground truth
#      exists, so only coverage and pairwise tool agreement are reported.
#
# Run from /data/projects/mpph_benchmark/ after run_all.sh's setup.
set -euo pipefail
cd "$(dirname "$0")/../.."   # /data/projects/mpph_benchmark

CPUS=28
KOFAM_DB=databases/kofam
EGGNOG_DB=databases/eggnog

source venv/bin/activate

annotate_all() {   # $1 = label, $2 = path to protein FASTA
  local label=$1 faa=$2
  echo "=== ${label}: mpph annotate ==="
  /usr/bin/time -v mpph annotate --fasta "$faa" --kofam-db "$KOFAM_DB" \
    --cpus "$CPUS" --out "results/mpph_${label}.tsv" \
    2> "results/mpph_${label}.time.log"

  echo "=== ${label}: KofamScan ==="
  /usr/bin/time -v singularity exec containers/kofamscan.sif exec_annotation \
    -p "${KOFAM_DB}/profiles" -k "${KOFAM_DB}/ko_list" -f mapper \
    -o "results/kofamscan_${label}.tsv" --cpu "$CPUS" "$faa" \
    2> "results/kofamscan_${label}.time.log"

  echo "=== ${label}: eggNOG-mapper ==="
  /usr/bin/time -v singularity exec containers/eggnog-mapper.sif emapper.py \
    -i "$faa" --itype proteins -m diamond --cpu "$CPUS" \
    --data_dir "$EGGNOG_DB" --output "emapper_${label}" --output_dir results/ \
    2> "results/emapper_${label}.time.log"
}

# --- 1. non-model KEGG genomes (ground truth available) --------------------
for org in vbs lbac; do
  annotate_all "$org" "genomes/${org}.faa"
done

# --- 2. real MAGs (no ground truth) ----------------------------------------
for mag in mags/*.faa; do
  annotate_all "$(basename "$mag" .faa)" "$mag"
done

# --- comparison -------------------------------------------------------------
echo "=== non-model genomes: accuracy vs KEGG ground truth ==="
for org in vbs lbac; do
  n_genes=$(grep -c "^>" "genomes/${org}.faa")
  python benchmarks/annotation/compare_annotation.py \
    --org "$org" --n-genes "$n_genes" \
    --mpph "results/mpph_${org}.tsv" \
    --kofamscan "results/kofamscan_${org}.tsv" \
    --emapper "results/emapper_${org}.emapper.annotations" \
    --ground-truth "ground_truth/link_ko_${org}.tsv" \
    --protein-id-map "ground_truth/protein_id_map_${org}.tsv" \
    --out "report/annotation_${org}.csv"
  echo
done

echo "=== MAGs: coverage + tool agreement only (no ground truth exists) ==="
for mag in mags/*.faa; do
  label=$(basename "$mag" .faa)
  n_genes=$(grep -c "^>" "$mag")
  python benchmarks/annotation/compare_annotation.py \
    --org "$label" --n-genes "$n_genes" \
    --mpph "results/mpph_${label}.tsv" \
    --kofamscan "results/kofamscan_${label}.tsv" \
    --emapper "results/emapper_${label}.emapper.annotations" \
    --out "report/annotation_${label}.csv"
  echo
done
