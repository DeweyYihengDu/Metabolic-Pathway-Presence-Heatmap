# Data sources and licensing

## KEGG

MPPH retrieves data at run time from the [KEGG REST API](https://www.kegg.jp/kegg/rest/keggapi.html)
(`list/genome`, `list/pathway`, `list/module`, `get/md:`, `link/ko`,
`get/br:br08901`, `info/kegg`, and `get/ko#####/kgml` for `mpph pathmap`'s
pathway/global-map diagrams). The KEGG release in effect for a run is recorded
in each run's `*_manifest.json`.

**The MIT License in this repository covers the MPPH source code only — it does
not cover KEGG data.** KEGG is a separate resource with its own terms:
academic-site use, service provision, and non-academic/commercial use are
treated differently. See the [KEGG copyright and disclaimer](https://www.kegg.jp/kegg/legal.html).

Practical implications:

- The local `.mpph_cache/` directory holds copies of KEGG responses for speed
  and reproducibility. **Do not redistribute it as a standalone KEGG mirror.**
- Before deploying MPPH as a public web service, verify that your use complies
  with the current KEGG terms.
- When you publish results, cite KEGG and record the KEGG release (the manifest
  captures it).

## KOfam

`mpph annotate` downloads KEGG's KOfam HMM profile database
(`ko_list.gz` + `profiles.tar.gz`, from KEGG's own distribution server) on
request, via `mpph annotate --setup-db DIR`. This is a separate KEGG-
distributed resource from the REST API above, with the same terms: **the
MIT License here covers MPPH's source code only, not KOfam's data.**

- The local `.mpph_kofam_db/` directory (default `--kofam-db`) holds the
  downloaded database. **Do not redistribute it as a standalone mirror**,
  the same caveat as `.mpph_cache/`.
- KO assignment via `mpph annotate` reimplements KofamScan's own scoring
  logic (Aramaki et al. 2020) against this database; if you publish results
  from it, cite KofamScan alongside MPPH.

## User-supplied annotations

`--user` input (KofamScan / eggNOG-mapper / DRAM KO tables, etc.) is your own
data; you are responsible for its licensing and for the accuracy of the
annotations you provide.

## Third-party data bundled as examples

`examples/Ecoli_*` (CPM-normalized expression, metadata, ranked list, GSEA output and
figures) is derived from a public, published dataset, included here under fair
academic use to demonstrate and validate `mpph gsea` against real data with a
known biological ground truth:

> Liou GG, Chao Kaberdina A, Wang WS, et al. "Combined Transcriptomic and
> Proteomic Profiling of *E. coli* under Microaerobic versus Aerobic
> Conditions..." *Int J Mol Sci.* 2022. doi:[10.3390/ijms23052570](https://doi.org/10.3390/ijms23052570).
> Data: [GEO GSE189154](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189154).

Gene identifiers were converted from the study's NCBI GeneIDs to KEGG KO ids via
KEGG's own `conv/eco/ncbi-geneid` and `link/ko/eco` endpoints (no third-party
mapping resource used). If you reuse these files, cite the original study above
in addition to MPPH.
