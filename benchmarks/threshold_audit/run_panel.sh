#!/usr/bin/env bash
# Sub-threshold collection across the multi-rank calibration panel.
#
# Runs ONE KOfam search over all panel proteomes concatenated, with gene ids
# prefixed by organism code, then splits the result back per organism. The
# alternative -- 90 separate runs -- would re-read the ~1.5 GB profile
# database 90 times for no benefit, since the search cost is otherwise
# linear in targets either way.
#
# `nice` throughout: this is a multi-hour job and is expected to share the
# machine with whatever else is running.
#
#   bash run_panel.sh panel/data databases/kofam 16
set -euo pipefail

DATA="${1:?usage: run_panel.sh PANEL_DATA_DIR KOFAM_DB CPUS}"
KOFAM="${2:?}"
CPUS="${3:-16}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

mkdir -p "$DATA/subthreshold"
COMBINED="$DATA/all_panel.faa"

if [[ ! -s "$COMBINED" ]]; then
  echo "[1/3] Concatenating proteomes with organism-prefixed ids ..."
  : > "$COMBINED"
  for faa in "$DATA"/genomes/*.faa; do
    org=$(basename "$faa" .faa)
    # `org|` prefix: the pipe never occurs in a RefSeq accession, so the
    # split back out afterwards is unambiguous.
    sed "s/^>/>${org}|/" "$faa" >> "$COMBINED"
  done
fi
echo "      $(grep -c '^>' "$COMBINED") protein(s) from $(ls "$DATA"/genomes/*.faa | wc -l) organism(s)"

if [[ ! -s "$DATA/subthreshold_all.tsv" ]]; then
  echo "[2/3] KOfam search (cpus=$CPUS, niced) -- expect several hours ..."
  nice -n 15 python "$ROOT/benchmarks/threshold_audit/collect_subthreshold.py" \
    --fasta "$COMBINED" --kofam-db "$KOFAM" --cpus "$CPUS" \
    --out "$DATA/subthreshold_all.tsv"
fi

echo "[3/3] Splitting per organism ..."
python - "$DATA" <<'PY'
import sys
from pathlib import Path
data = Path(sys.argv[1])
out_dir = data / "subthreshold"
out_dir.mkdir(exist_ok=True)
header = "gene\tko\tscore\tthreshold\tdelta\tscore_type\n"
handles = {}
try:
    with open(data / "subthreshold_all.tsv", encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            org, rest = line.split("|", 1)
            if org not in handles:
                handles[org] = open(out_dir / f"subthreshold_{org}.tsv", "w",
                                    encoding="utf-8")
                handles[org].write(header)
            handles[org].write(rest)
finally:
    for h in handles.values():
        h.close()
print(f"wrote {len(handles)} per-organism file(s) to {out_dir}")
PY

echo "=== PANEL RUN COMPLETE ==="
