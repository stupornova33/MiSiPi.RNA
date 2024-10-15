#' function to calculate GC content of reads
#' takes two data tables of reads with sequences
#' returns GC percent
#'
#' @param seqs a vector of sequences
#' @return GC percent

#' @export

get_GC_content <- function(seqs){
   all_dat <- Biostrings::DNAStringSet(seqs)
   freqs <- Biostrings::alphabetFrequency(all_dat, baseOnly = T, as.prob = T)
   GC <- (sum(freqs[,2]) + sum(freqs[,3]))/nrow(freqs)
   all_dat <- NULL
   freqs <- NULL
   return(GC)
}
