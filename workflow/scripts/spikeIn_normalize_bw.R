# setup --------------------------------------------------------------------------------------------------------
library(GenomicRanges)
library(rtracklayer)
# library(BSgenome.Dmelanogaster.UCSC.dm6)
library(tidyverse)

# import data --------------------------------------------------------------------------------------------------
# get filename of bigwig that will be used for normalization
bw <- snakemake@input[["bw"]]

# read in normalization factors for merged replicates
scaling_factors <- read_tsv(snakemake@input[["scaling_factors"]])

# get sample name to use for finding the correct scaling factor
file_bn <- gsub(".bw", "", basename(bw))

# rescale bigwig file based on spike-in normalization factor ---------------------------------------------------

# get scaling factor
scaling_factor <- scaling_factors %>%
  filter(IP_group == file_bn) %>%
  select(norm_scaling_factor) %>%
  pull()


# normalize bigwig
scaled.gr <- import(bw) %>%
  as.data.frame() %>%
  mutate(score = score * scaling_factor) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


# seqinfo <- seqinfo(BSgenome.Dmelanogaster.UCSC.dm6)
# keepers <- seqlevels(BSgenome.Dmelanogaster.UCSC.dm6)[seqlevels(BSgenome.Dmelanogaster.UCSC.dm6) %in% seqlevels(scaled.gr)]
# seqinfo <- keepSeqlevels(seqinfo, keepers)
# seqinfo(scaled.gr) <- seqinfo

# export normalized bigwig -----------------------------------------------------
export(scaled.gr, snakemake@output[["bw"]], format = "bw")


