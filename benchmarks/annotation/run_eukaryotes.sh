#!/usr/bin/env bash
# Tier 4: eukaryotes. Same three tools, same KOfam database, real KEGG
# ground truth. Yeast is single-isoform (6,021 proteins == 6,021 KEGG
# mappings); Arabidopsis is not (48,265 vs 27,562), which is what the
# evaluable-set restriction in compare_annotation.py exists for.
set -euo pipefail
cd "$(dirname "$0")/../.."   # the benchmark root, alongside run_all.sh

CPUS=28; KOFAM_DB=databases/kofam; EGGNOG_DB=databases/eggnog
source venv/bin/activate

for org in sce ath; do
  faa="genomes/${org}.faa"
  echo "=== ${org}: mpph annotate ($(grep -c '^>' $faa) proteins) ==="
  /usr/bin/time -v mpph annotate --fasta "$faa" --kofam-db "$KOFAM_DB" \
    --cpus "$CPUS" --out "results/mpph_${org}.tsv" 2> "results/mpph_${org}.time.log"
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

for org in sce ath; do
  n=$(grep -c "^>" "genomes/${org}.faa")
  python benchmarks/annotation/compare_annotation.py --org "$org" --n-genes "$n" \
    --mpph "results/mpph_${org}.tsv" --kofamscan "results/kofamscan_${org}.tsv" \
    --emapper "results/emapper_${org}.emapper.annotations" \
    --ground-truth "ground_truth/link_ko_${org}.tsv" \
    --protein-id-map "ground_truth/protein_id_map_${org}.tsv" \
    --out "report/annotation_${org}.csv"
  echo
done
