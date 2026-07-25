# Input formats

## KEGG sources

- **Taxon name** (positional): `mpph Vibrio`. Matching is by organism name:
  `--match word` (default, whole word), `prefix`, `exact`, or `regex`. Reliable
  for genus/species, not for family/phylum. Preview with `--dry-run`.
- **Organism codes**: `--codes codes.txt` (one KEGG code per line, `#` comments).

## Your own genomes / MAGs (`--user`)

`--user` accepts a directory (one file per genome, named by its stem; a
duplicate stem across extensions, e.g. `sample.txt` and `sample.tsv`, is
rejected rather than silently overwritten) or a single file. `--input-format`:

- `auto` (default) — scan any `K#####` id in the whole file (one sample).
  For a two-column file, also detects a genuine `sample<TAB>KO` long table --
  but only when a sample id actually repeats across rows. A file with one
  unique id per row (e.g. KofamScan `--format mapper`, `gene<TAB>KO`) is
  indistinguishable from a one-sample table by structure alone, so it is read
  as a single sample instead of silently turning every gene into its own
  "sample".
- `ko-list` — force KO-id scanning of the whole file as one sample (KofamScan
  `--format mapper` output works).
- `long` — force `sample<TAB>KO` long-table parsing, for the rare genuine
  multi-sample table where every sample happens to contribute exactly one row
  (the one case `auto` cannot infer on its own).
- `eggnog` — read the `KEGG_ko` column of an eggNOG-mapper `.annotations` file.

Don't have KO annotations yet? `mpph annotate --fasta genome.faa` (see
[Subcommands](subcommands.md#mpph-annotate)) produces this exact
`gene<TAB>K#####` mapper shape directly from a protein FASTA -- no external
KofamScan/eggNOG-mapper install needed. Point `--user` at its output file
directly, or at a directory of several genomes' worth of `mpph annotate`
output (one file per genome) for the standard one-file-per-sample layout.

## Sample sheets (reproducible)

A TSV with `sample_id` and `annotation_file`, optional `input_format`, and any
metadata columns:

```
sample_id   annotation_file        input_format   group     completeness
MAG001      data/MAG001.kofam.tsv  ko-list        surface   95.4
MAG002      data/MAG002.emapper    eggnog         deep      82.3
```

Validate with `mpph validate --samples samples.tsv`. Metadata columns feed
`compare` / `ordination` grouping. `sample_id` and `annotation_file` must be
non-blank and `sample_id` must be unique -- a blank `annotation_file` is
rejected rather than silently resolving to the sheet's base directory itself.

## Genome QC metadata

`mpph.samplesheet.import_checkm2()` and `import_gtdbtk()` read CheckM2
`quality_report.tsv` (completeness/contamination) and GTDB-Tk `summary.tsv`
(taxonomy), keyed by genome, ready to join onto a sample sheet.
