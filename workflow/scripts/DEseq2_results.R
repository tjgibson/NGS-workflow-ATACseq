# setup ------------------------------------------------------------------------
## load necessary packages
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))

# import dds object ------------------------------------------------------------
## import dds2 object
dds <- readRDS(snakemake@input[["dds"]])

# get DEseq results for all contrasts ------------------------------------------
contrast <- c("condition", snakemake@params[["contrast"]])
results <- results(dds, contrast = contrast, alpha = snakemake@params[["padj_cutoff"]]) %>% 
  as.data.frame() %>% 
  arrange(padj) %>% 
  rownames_to_column(var = "peak_id")

# add additional information to results table ----------------------------------
# #read in peak file
peaks <- rtracklayer::import(snakemake@input[["peaks"]]) %>% 
  as.data.frame() %>%
  dplyr::select(seqnames, start, end, name, signalValue, pValue, qValue) %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end, peak_id = name, MACS2_enrichment = signalValue, MACS2_pValue = pValue, MACS2_qValue = qValue)

# add gene symbol to results
out_table <- dplyr::left_join(results, peaks, by = "peak_id")

## add column indicating if gene is differentially expressed with padj < 0.05 and FC > 2
out_table <- out_table %>%
  dplyr::mutate(is_diff = (padj < snakemake@params[["padj_cutoff"]] & (abs(log2FoldChange) > snakemake@params[["FC_cutoff"]]))) %>%
  replace_na(list(is_diff = FALSE))

# write output file ------------------------------------------------------------
write_tsv(out_table, snakemake@output[[1]])