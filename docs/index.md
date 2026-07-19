# MPPH documentation

MPPH turns KEGG (or your own KO annotations) into auditable, comparable
microbial metabolic-function profiles.

- [Methods](methods.md) — what "presence" and "completeness" mean, how modules
  are scored, and how the dendrogram should (and should not) be interpreted.
- [Input formats](input-formats.md) — taxon / codes / user KO annotations,
  eggNOG, sample sheets, CheckM2 / GTDB-Tk.
- [Subcommands](subcommands.md) — `run`, `traits`, `explain`, `compare`, `pan`,
  `ordination`, `community`, `report`, `validate`, `enrich`.
- [Output schema](outputs.md) — every file MPPH writes.

Quick start:

```bash
pip install .
mpph Prochlorococcus --cluster --drop-core          # pathway-map heatmap
mpph Prochlorococcus --completeness --cluster        # module completeness
mpph traits Prochlorococcus --panel biogeochemistry  # trait panel
mpph enrich --study genes.txt --background-organism pmt \
            --ontology kegg-pathway                  # enrichment (ORA)
```

See the top-level `README.md` for the gallery and the public Python API.
