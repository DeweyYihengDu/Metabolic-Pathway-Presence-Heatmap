#!/usr/bin/env bash
# Fetch proteomes + KEGG ground truth for the multi-rank calibration panel.
#
# The seven-genome benchmark cannot answer whether calibrating on a closer
# relative helps -- it has exactly 2 same-genus pairs. This panel has 144.
#
# Two inputs per organism:
#   * proteome FASTA   -- from NCBI RefSeq, resolved via assembly_summary
#   * link/ko/<org>    -- KEGG's own curated KO assignments (the ground truth)
#   * conv/<org>/ncbi-proteinid -- id namespace bridge between the two
#
# Organisms with no RefSeq protein FASTA or no conv table are skipped and
# listed at the end rather than silently dropped: a panel that quietly loses
# its close pairs would answer the wrong question.
#
#   bash fetch_panel.sh panel.tsv panel_data/
set -uo pipefail

PANEL="${1:?usage: fetch_panel.sh PANEL.tsv OUTDIR}"
OUT="${2:?usage: fetch_panel.sh PANEL.tsv OUTDIR}"
mkdir -p "$OUT"/{genomes,ground_truth,logs}

SUMMARY="$OUT/assembly_summary_refseq.txt"
if [[ ! -s "$SUMMARY" ]]; then
  echo "[1/3] Fetching NCBI RefSeq assembly summary (~250 MB) ..."
  wget -q -O "$SUMMARY" \
    https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
fi

skipped=()
n_ok=0
total=$(( $(wc -l < "$PANEL") - 1 ))
i=0

echo "[2/3] Fetching $total organism(s) ..."
while IFS=$'\t' read -r org name kingdom domain clade genus; do
  [[ "$org" == "org" ]] && continue
  i=$((i + 1))
  faa="$OUT/genomes/${org}.faa"
  gt="$OUT/ground_truth/link_ko_${org}.tsv"
  cv="$OUT/ground_truth/protein_id_map_${org}.tsv"

  # KEGG side first: an organism with no curated KO set is useless here, and
  # finding that out before a 100 MB download saves the bandwidth.
  [[ -s "$gt" ]] || curl -sf "https://rest.kegg.jp/link/ko/${org}" -o "$gt"
  [[ -s "$cv" ]] || curl -sf "https://rest.kegg.jp/conv/${org}/ncbi-proteinid" -o "$cv"
  if [[ ! -s "$gt" || ! -s "$cv" ]]; then
    skipped+=("$org (no KEGG ground truth or conv table)")
    continue
  fi

  if [[ ! -s "$faa" ]]; then
    # KEGG records the source assembly; resolve it to an FTP path via the
    # summary rather than guessing the directory layout.
    acc=$(curl -sf "https://rest.kegg.jp/get/gn:${org}" \
          | awk '/DATA_SOURCE/ {for(j=1;j<=NF;j++) if($j ~ /^GC[AF]_/) {print $j; exit}}')
    if [[ -z "${acc:-}" ]]; then
      skipped+=("$org (no assembly accession in KEGG record)"); continue
    fi
    # KEGG's DATA_SOURCE usually names the *GenBank* assembly (GCA_...), while
    # column 1 of the RefSeq summary holds GCF_ accessions -- the GCA pairing
    # lives in column 18. Matching only column 1 skipped ~90% of the panel.
    # Compared version-stripped and exactly, not by regex: a substring match
    # on an accession can hit the wrong assembly.
    ftp=$(awk -F'\t' -v a="${acc%%.*}" '
        {split($1, x, "."); split($18, y, ".");
         if (x[1] == a || y[1] == a) {print $20; exit}}' "$SUMMARY")
    if [[ -z "${ftp:-}" ]]; then
      skipped+=("$org ($acc not in RefSeq summary)"); continue
    fi
    base=$(basename "$ftp")
    if ! wget -q -O "${faa}.gz" "${ftp}/${base}_protein.faa.gz"; then
      rm -f "${faa}.gz"; skipped+=("$org (no protein.faa.gz)"); continue
    fi
    gunzip -f "${faa}.gz" || { skipped+=("$org (gunzip failed)"); continue; }
  fi
  n_ok=$((n_ok + 1))
  printf '\r  %d/%d ok=%d skipped=%d' "$i" "$total" "$n_ok" "${#skipped[@]}"
done < "$PANEL"
echo

echo "[3/3] $n_ok organism(s) ready in $OUT/genomes"
if ((${#skipped[@]})); then
  printf '%s\n' "${skipped[@]}" > "$OUT/logs/skipped.txt"
  echo "${#skipped[@]} skipped -- see $OUT/logs/skipped.txt"
fi
