# function to extract all reads in specified region on a given strand
# takes bam object, chrom_name, reg_start, and reg_stop
# extracts rname, pos, qwidth, and seq from bam obj
#
# @param bam_obj a bam obj created by open_bam
# @param chrom_name a name passed in by user/bed file
# @param reg_start a whole number
# @param reg_stop a whole number
# @param strand - Character - "plus" or "minus"
# @return chrom obj

.get_chr <- function(bam_obj, chrom_name, reg_start, reg_stop, strand = c("minus", "plus", "-", "+")) {
  strand <- match.arg(strand)
  is_minus <- switch(strand,
    "minus" = TRUE,
    "plus" = FALSE,
    "-" = TRUE,
    "+" = FALSE
  )

  Rsamtools::open.BamFile(bam_obj)
  which <- GenomicRanges::GRanges(seqnames = chrom_name, IRanges::IRanges(reg_start, reg_stop))
  bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = is_minus, isUnmappedQuery = FALSE), what = c("rname", "pos", "qwidth", "seq"), which = which)
  result <- Rsamtools::scanBam(bam_obj, param = bam_scan)[[1]]
  Rsamtools::close.BamFile(bam_obj)
  return(result)
}
