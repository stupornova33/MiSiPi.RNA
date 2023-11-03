#' create an object from read data from plus strand
#' takes bam obj, chrom name, reg_start, and reg_stop
#'
#' @param bam_obj a bam object created by openBamFile
#' @param reg_start an integer
#' @param reg_stop an integer
#' @return max_overhang



#' @export
getChrPlus <- function(bam_obj, chrom_name, reg_start, reg_stop){
   Rsamtools::open.BamFile(bam_obj)
   which <- GenomicRanges::GRanges(seqnames=c(chrom_name), IRanges::IRanges(c(reg_start), c(reg_stop)))
   bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = FALSE, isUnmappedQuery = FALSE), what=c('rname', 'pos', 'qwidth', 'seq'), which=which)
   result <- Rsamtools::scanBam(bam_obj, param=bam_scan)[[1]]
   Rsamtools::close.BamFile(bam_obj)
   return(result)
}
