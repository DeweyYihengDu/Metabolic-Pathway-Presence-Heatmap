#!/usr/bin/env bash
# Sub-threshold precision audit across every genome with KEGG ground truth.
set -euo pipefail
cd "$(dirname "$0")/../.."
source venv/bin/activate
CPUS=24
for org in bsu mja vbs lbac sce ath; do
  echo "=== ${org}: collecting ==="
  /usr/bin/time -v python benchmarks/threshold_audit/collect_subthreshold.py \
    --fasta "genomes/${org}.faa" --kofam-db databases/kofam --cpus "$CPUS" \
    --out "threshold_audit/subthreshold_${org}.tsv" \
    2> "threshold_audit/collect_${org}.time.log"
  echo "=== ${org}: scoring ==="
  python benchmarks/threshold_audit/score_precision.py \
    --hits "threshold_audit/subthreshold_${org}.tsv" --org "$org" \
    --ground-truth "ground_truth/link_ko_${org}.tsv" \
    --protein-id-map "ground_truth/protein_id_map_${org}.tsv" \
    --out "threshold_audit/precision_${org}.csv"
done
echo "=== AUDIT COMPLETE ==="
