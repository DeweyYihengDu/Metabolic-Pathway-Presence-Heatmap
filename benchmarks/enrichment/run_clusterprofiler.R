#!/usr/bin/env Rscript
# Runs clusterProfiler's GSEA() and enricher() against the *same* ranked
# list / TERM2GENE category definitions mpph gsea/enrich used, so the
# comparison isolates the statistical implementation, not the database.
# Intended to run inside the bioconductor-clusterprofiler container:
#   singularity exec containers/clusterprofiler.sif Rscript \
#     benchmarks/enrichment/run_clusterprofiler.R
suppressMessages(library(clusterProfiler))

args <- commandArgs(trailingOnly = TRUE)
ranked_path   <- if (length(args) >= 1) args[1] else "examples/Ecoli_ranked_logFC.tsv"
term2gene_path <- if (length(args) >= 2) args[2] else "term2gene_metabolism.tsv"
term2name_path <- if (length(args) >= 3) args[3] else "term2name.tsv"
out_dir       <- if (length(args) >= 4) args[4] else "report"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ranked_df <- read.delim(ranked_path, header = FALSE,
                        col.names = c("gene", "score"))
gene_list <- setNames(ranked_df$score, ranked_df$gene)
gene_list <- sort(gene_list, decreasing = TRUE)  # GSEA requires descending order

# colClasses = "character" is load-bearing: KEGG pathway ids are 5-digit
# zero-padded strings (e.g. "01210"), and read.delim()'s own type inference
# reads an all-numeric-looking column as integer by default, silently
# stripping the leading zero ("01210" -> 1210 -> "1210") -- the exact same
# class of bug pandas' auto dtype inference causes on this same id shape
# (see compare_annotation.py's/mpph's own report.py history). Confirmed by
# running without this once: clusterProfiler's own GSEA() output ID column
# came back as "1210" for what should have been "01210", so it never matched
# mpph's category_id and the whole comparison silently found zero overlap.
term2gene <- read.delim(term2gene_path, header = FALSE,
                        colClasses = "character", col.names = c("term", "gene"))
term2name <- read.delim(term2name_path, header = FALSE,
                        colClasses = "character", col.names = c("term", "name"))

# pvalueCutoff = 1: report every tested pathway, not just significant ones --
# matching mpph gsea's own behaviour (it reports the whole tested universe,
# with q_value as a column to filter on afterward, not a hard cutoff).
# minGSSize/maxGSSize = 15/500: clusterProfiler's own defaults are 10/500,
# NOT mpph gsea's (--min-size 15 --max-size 500) -- left at clusterProfiler's
# default the first time this ran, 68 pathways got tested here vs mpph's 55,
# entirely from the size-filter mismatch, not a real GSEA-implementation
# difference. Matched explicitly so both tools test the identical pathway
# set and the comparison isolates the statistic, not an accidental parameter
# mismatch.
gsea_result <- GSEA(gene_list, TERM2GENE = term2gene, TERM2NAME = term2name,
                    minGSSize = 15, maxGSSize = 500,
                    pvalueCutoff = 1, seed = TRUE, verbose = FALSE)
gsea_df <- as.data.frame(gsea_result)
write.csv(gsea_df, file.path(out_dir, "clusterprofiler_gsea.csv"), row.names = FALSE)
cat(sprintf("GSEA: %d pathways tested. Wrote %s\n", nrow(gsea_df),
           file.path(out_dir, "clusterprofiler_gsea.csv")))

# ORA: derive a study set the same simple way any DE threshold would --
# |log2FC| in the top/bottom 10% of the ranking, background = every ranked gene.
# (This mirrors going from a continuous DE result to a discrete "significant
# gene" list; it does not require any additional input data beyond what
# mpph gsea already used.)
cutoff_hi <- quantile(gene_list, 0.90)
cutoff_lo <- quantile(gene_list, 0.10)
study <- names(gene_list)[gene_list >= cutoff_hi | gene_list <= cutoff_lo]
universe <- names(gene_list)

ora_result <- enricher(study, universe = universe, TERM2GENE = term2gene,
                       TERM2NAME = term2name, pvalueCutoff = 1, minGSSize = 2)
ora_df <- as.data.frame(ora_result)
write.csv(ora_df, file.path(out_dir, "clusterprofiler_ora.csv"), row.names = FALSE)
cat(sprintf("ORA: %d study genes, %d pathways tested. Wrote %s\n",
           length(study), nrow(ora_df),
           file.path(out_dir, "clusterprofiler_ora.csv")))

# Also export the study/universe sets so mpph enrich is tested on the exact
# same gene lists, not a different DE threshold.
writeLines(study, file.path(out_dir, "study_genes.txt"))
writeLines(universe, file.path(out_dir, "background_genes.txt"))
