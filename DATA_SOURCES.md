# Data sources and licensing

## KEGG

MPPH retrieves data at run time from the [KEGG REST API](https://www.kegg.jp/kegg/rest/keggapi.html)
(`list/genome`, `list/pathway`, `list/module`, `get/md:`, `link/ko`,
`get/br:br08901`, `info/kegg`). The KEGG release in effect for a run is recorded
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

## User-supplied annotations

`--user` input (KofamScan / eggNOG-mapper / DRAM KO tables, etc.) is your own
data; you are responsible for its licensing and for the accuracy of the
annotations you provide.
