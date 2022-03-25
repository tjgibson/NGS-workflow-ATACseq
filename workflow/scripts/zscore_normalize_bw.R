# setup -----------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))
              

zscore_bw <- function(bw) {
  require(tidyverse)
  require(rtracklayer)
  require(GenomicRanges)
  
  # import bigwig file to Granges
  if (typeof(bw) == "character") {
    message("reading bigwig file")
    bw <- import(bw)
  }
  
  # if using a spike-in, filter the seqlevels to only the reference genome
  if (snakemake@config[["use_spikeIn"]]) {
    message("removing spikeIn chromosomes")
    ref_chroms <- seqlevels(bw)[!grepl("spikeIn_", seqlevels(bw))]
    bw <- keepSeqlevels(bw, ref_chroms, pruning.mode = "coarse")
  }
  
  if (snakemake@config[["filter_chroms"]]) {
    message("filtering reference chromosomes")
    keep_chroms <- read_tsv(snakemake@config[["keep_chroms"]], col_names = c("chromosome"))
    ref_chroms <- seqlevels(bw)[seqlevels(bw) %in% keep_chroms$chromosome]
    bw <- keepSeqlevels(bw, ref_chroms, pruning.mode = "coarse")
  }
  
  
  # for large regions with the same score, expand into equal sized bins
  message("binning genome")
  min_binsize <- min(width(bw))
  all_bins <- tileGenome(seqinfo(bw), tilewidth=min_binsize,cut.last.tile.in.chrom=TRUE)
  
  message("getting scores for all bins")
  # add the coverage/score for both input and IP
  all_bins <- subsetByOverlaps(all_bins, bw)
  overlaps <- findOverlaps(all_bins, bw)
  all_bins$score[overlaps@from] <- bw$score[overlaps@to]
  
  # perform z-score normalization
  message("performing z-score normalization")
  all_bins$zscore <- scale(all_bins$score)[,1]
  all_bins$score <- NULL
  all_bins$score <- all_bins$zscore
  all_bins$zscore <- NULL
  # collapse adjacent bins with same score
  collapsed <- unlist(GenomicRanges::reduce(split(all_bins, ~score)))
  collapsed$score <- as.numeric(names(collapsed))
  names(collapsed) <- NULL
  all_bins <- collapsed
  
  #set seqinfo for z-score normalized version
  seqinfo(all_bins) <- seqinfo(bw)
  
  return(all_bins)
}


# perform z-score normalization and write new bigwig files ---------------------
zscore.gr <- zscore_bw(snakemake@input[[1]])
export(zscore.gr, snakemake@output[[1]])