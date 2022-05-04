# setup -----------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))

# define functions -------------------------------------------------------------
# function to extract peak summits from narrowPeak file
# peak start and end values will be replaced with the summit location
extract_summits <- function(gr, extend_width = 1L) {
  # verify input is a GRanges object
  if (!inherits(gr, "GRanges") ) {
    stop("x must be a GRanges object")
  }
  # verify that gr object is in narrowPeak format
  if ( !all(names(mcols(gr))  %in% c("name", "score", "signalValue", "pValue",  "qValue",  "peak"))) {
    stop(strwrap("GRanges object does not appear to be in narrowPeak format. Object should contain the following metadata columns: name, score, signalValue, pValue, qValue and peak"))
  }
  
  if (all(gr$peak == -1)) {
    stop("All values for 'peak' column == -1, indicating that summits were not called. Verify that your peak caller called peak summits")
  }
  
  # replace peak start and end values with summit location
  summit <- start(gr) + gr$peak - 1
  start(gr) <- summit
  end(gr) <- summit
  
  # resize peaks to desired width
  gr <- resize(gr, width = extend_width, fix = "center")

  # remove now meaningless peak column
  gr$peak <- NULL
  
  return(gr)
}


# read peaks -------------------------------------------------------------------
peak_fn <- snakemake@input[[1]]
peaks <- rtracklayer::import(peak_fn)


# extend peak summits ----------------------------------------------------------
extended_summits <- peaks %>% 
  extract_summits(extend_width = as.integer(snakemake@params[["extend_width"]]))

# export filtered peaks to file  -----------------------------------------------
extended_summits %>% 
  as.data.frame() %>% 
  select(1:3, 6:7, 5, 8:10) %>% 
  mutate(strand = ".") %>% 
  mutate(peak = -1) %>% 
  write_tsv(snakemake@output[[2]], col_names = FALSE)

rtracklayer::export(extended_summits, snakemake@output[[1]])

# export SAF file for featureCounts --------------------------------------------
out_peaks.df <- as.data.frame(extended_summits)
out_peaks.df$strand <- "."
peaks.saf <- out_peaks.df %>% 
  select(c("name", "seqnames", "start", "end", "strand"))
names(peaks.saf) <- c( "GeneID", "Chr", "Start", "End", "Strand")

write_tsv(peaks.saf, snakemake@output[[3]])
