# setup ------------------------------------------------------------------------
# load required packages
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

# import data ------------------------------------------------------------------
# read in all relevant stat files
mapping_stats <- snakemake@input[["mapping_stats"]] %>%
  map_dfr(read_tsv, .id = "source", col_names = c("chromosome", "chromosome_size", "mapped_reads", "unmapped_reads")) %>%
  mutate(sample = basename(source))

# add column marking reference or spike-in chromosomes
mapping_stats <- mapping_stats %>%
  mutate(reference_genome = case_when(str_detect(chromosome, "spikeIn")  ~ "spikeIn",
                                      !str_detect(chromosome, "spikeIn")  ~ "reference")) 



mapping_stats <- mapping_stats %>% 
  group_by(sample, reference_genome) %>% 
  summarise(total_reads = sum(mapped_reads)) %>%
  pivot_wider(names_from = reference_genome, values_from = total_reads) %>%
  mutate(total_mapped_reads = reference + spikeIn) %>%
  mutate(percent_reference = reference / total_mapped_reads * 100, percent_spikeIn = spikeIn / total_mapped_reads * 100)


# compute normalization factors ------------------------------------------------


# rescale bigwig file based on spike-in normalization factor -------------------
