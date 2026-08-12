#!/usr/bin/env bash
# Runs mpph annotate, KofamScan, and eggNOG-mapper on the same protein FASTA
# for each of eco/bsu/mja, then compare_annotation.py on the results.
# Intended to run from /data/projects/mpph_benchmark/ on Hermes, after
# setup_containers.log and setup_mpph_kofam.log both show their DONE markers,
# and after the eggNOG-mapper database has been downloaded into databases/eggnog.
set -euo pipefail
cd "$(dirname "$0")/../.."   # /data/projects/mpph_benchmark

CPUS=28   # leave a few of Hermes' 32 cores free for other work
KOFAM_DB=databases/kofam
EGGNOG_DB=databases/eggnog

source venv/bin/activate

for org in eco bsu mja; do
  faa="genomes/${org}.faa"
  echo "=== ${org}: mpph annotate ==="
  /usr/bin/time -v mpph annotate --fasta "$faa" --kofam-db "$KOFAM_DB" \
    --cpus "$CPUS" --out "results/mpph_${org}.tsv" \
    2> "results/mpph_${org}.time.log"

  echo "=== ${org}: KofamScan ==="
  /usr/bin/time -v singularity exec containers/kofamscan.sif exec_annotation \
    -p "${KOFAM_DB}/profiles" -k "${KOFAM_DB}/ko_list" -f mapper \
    -o "results/kofamscan_${org}.tsv" --cpu "$CPUS" "$faa" \
    2> "results/kofamscan_${org}.time.log"

  echo "=== ${org}: eggNOG-mapper ==="
  /usr/bin/time -v singularity exec containers/eggnog-mapper.sif emapper.py \
    -i "$faa" --itype proteins -m diamond --cpu "$CPUS" \
    --data_dir "$EGGNOG_DB" --output "emapper_${org}" --output_dir results/ \
    2> "results/emapper_${org}.time.log"
done

echo "=== spot-check known genes before trusting full results ==="
for org in eco bsu mja; do
  echo "--- $org: first 5 mpph calls ---"
  head -5 "results/mpph_${org}.tsv"
done

echo "=== comparison ==="
for org in eco bsu mja; do
  n_genes=$(grep -c "^>" "genomes/${org}.faa")
  python benchmarks/annotation/compare_annotation.py \
    --org "$org" --n-genes "$n_genes" \
    --mpph "results/mpph_${org}.tsv" \
    --kofamscan "results/kofamscan_${org}.tsv" \
    --emapper "results/emapper_${org}.emapper.annotations" \
    --ground-truth "ground_truth/link_ko_${org}.tsv" \
    --protein-id-map "ground_truth/protein_id_map_${org}.tsv" \
    --out "report/annotation_${org}.csv"
done
