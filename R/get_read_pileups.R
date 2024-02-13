#' A function to extract read pileups at specified start and stop
#' @param start a whole number
#' @param stop a whole number
#' @param bam_scan a bam scan object
#' @param input_file a string

#' @return a data table
#' @export

get_read_pileups <- function(start, stop, bam_scan, input_file){
   seqnames <- pos <- count <- NULL
   params <- Rsamtools::PileupParam(max_depth=10000, min_base_quality=0, min_mapq=0, min_nucleotide_depth=0, distinguish_strands=TRUE,
                                    distinguish_nucleotides=TRUE, ignore_query_Ns=TRUE, include_deletions=FALSE, include_insertions=FALSE, left_bins=NULL,
                                    query_bins=NULL, cycle_bins=NULL)

   pileups <- Rsamtools::pileup(input_file, index=(stringr::str_c(input_file, '','.bai')), scanBamParam=bam_scan, pileupParam=params) %>%
      dplyr::select(-c(seqnames, strand))

   dt <- pileups %>% dplyr::group_by(pos) %>% dplyr::summarise(count = sum(count))
   return(dt)
}
