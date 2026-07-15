# Input formats

## KEGG sources

- **Taxon name** (positional): `mpph Vibrio`. Matching is by organism name:
  `--match word` (default, whole word), `prefix`, `exact`, or `regex`. Reliable
  for genus/species, not for family/phylum. Preview with `--dry-run`.
- **Organism codes**: `--codes codes.txt` (one KEGG code per line, `#` comments).

## Your own genomes / MAGs (`--user`)

`--user` accepts a directory (one file per genome, named by its stem) or a
single file. `--input-format`:

- `auto` (default) — scan any `K#####` id; also detects a `sample<TAB>KO`
  long table.
- `ko-list` — force KO-id scanning (KofamScan `--format mapper` output works).
- `eggnog` — read the `KEGG_ko` column of an eggNOG-mapper `.annotations` file.

## Sample sheets (reproducible)

A TSV with `sample_id` and `annotation_file`, optional `input_format`, and any
metadata columns:

```
sample_id   annotation_file        input_format   group     completeness
MAG001      data/MAG001.kofam.tsv  ko-list        surface   95.4
MAG002      data/MAG002.emapper    eggnog         deep      82.3
```

Validate with `mpph validate --samples samples.tsv`. Metadata columns feed
`compare` / `ordination` grouping.

## Genome QC metadata

`mpph.samplesheet.import_checkm2()` and `import_gtdbtk()` read CheckM2
`quality_report.tsv` (completeness/contamination) and GTDB-Tk `summary.tsv`
(taxonomy), keyed by genome, ready to join onto a sample sheet.
