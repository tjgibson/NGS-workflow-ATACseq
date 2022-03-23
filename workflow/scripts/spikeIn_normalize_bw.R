# setup --------------------------------------------------------------------------------------------------------
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)

# import data --------------------------------------------------------------------------------------------------
# get filename of bigwig that will be used for normalization
bw <- snakemake@input[["bw"]]

# read in normalization factors for merged replicates
scaling_factors <- read_tsv(snakemake@input[["scaling_factors"]])

# get sample name to use for finding the correct scaling factor
file_bn <- gsub(".bw", "", basename(bw))

## for testing only ##
print("bw variable:")
type(bw)
str(bw)
print(bw)
print(" ")

print("scaling_factors variable:")
type(scaling_factors)
str(scaling_factors)
print(scaling_factors)
print(" ")
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

# set seqinfo of normalized bigwigs
seqlevels(scaled.gr) <- seqlevels(BigWigFile(bw))
seqinfo(scaled.gr) <- seqinfo(BigWigFile(bw))


# export normalized bigwig -----------------------------------------------------
export(scaled.gr, snakemake@output[["bw"]], format = "bw")


