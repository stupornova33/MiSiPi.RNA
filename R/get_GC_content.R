#' function to calculate GC content of reads
#' takes two data tables of reads with sequences
#' returns GC percent
#'
#' @param dt1 a data table of read data
#' @param dt2 a data table of read data
#' @return GC percent



#' @export


get_GC_content <- function(data){
   all_dat <- Biostrings::DNAStringSet(data$seq)
   freqs <- Biostrings::alphabetFrequency(all_dat, baseOnly = T, as.prob = T)
   GC <- (sum(freqs[,2]) + sum(freqs[,3]))/nrow(freqs)
   return(GC)
}
