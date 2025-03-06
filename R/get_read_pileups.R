# A function to extract read pileups at specified start and stop
# @param chrom_name chromosome name
# @param start a whole number
# @param stop a whole number
# @param strand either + or -
# @param bam_file path to bam file
# @return a data frame

.get_read_pileups <- function(chrom_name, start, stop, strand, bam_file) {
  seqnames <- pos <- count <- NULL
  
  is_minus <- switch(
    strand,
    "minus" = TRUE,
    "plus" = FALSE,
    "-" = TRUE,
    "+" = FALSE
  )
  
  which <- GenomicRanges::GRanges(
    seqnames = chrom_name,
    IRanges::IRanges(start, stop)
  )
  
  scan_params <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(isMinusStrand = is_minus),
    what = c("rname", "pos", "qwidth"),
    which = which
  )
  
  pileup_params <- Rsamtools::PileupParam(
    max_depth = 10000,
    min_base_quality = 0,
    min_mapq = 0,
    min_nucleotide_depth = 0,
    distinguish_strands = TRUE,
    distinguish_nucleotides = TRUE,
    ignore_query_Ns = TRUE,
    include_deletions = FALSE,
    include_insertions = FALSE,
    left_bins = NULL,
    query_bins = NULL,
    cycle_bins = NULL
  )

  pileups <- Rsamtools::pileup(
    bam_file,
    index = (stringr::str_c(bam_file, "", ".bai")),
    scanBamParam = scan_params,
    pileupParam = pileup_params) %>%
    dplyr::select(-c(seqnames, strand))

  dt <- pileups %>%
    dplyr::group_by(pos) %>%
    dplyr::summarise(count = sum(count))
  
  return(dt)
}
