#' function to extract all reads in specified region on minus strand
#' takes bam object, chrom_name, reg_start, and reg_stop
#' extracts rname, pos, qwidth, and seq from bam obj
#' 
#' 
#' @param bam_obj a bam obj created by openBamFile
#' @param chrom_name a name passed in by user/bed file
#' @param reg_start a whole number
#' @param reg_stop a whole number
#' @return chrom obj



#' @export

getChrMinus <- function(bam_obj, chrom_name, reg_start, reg_stop){
   Rsamtools::open.BamFile(bam_obj)
   which <- GenomicRanges::GRanges(seqnames=chrom_name, IRanges::IRanges(reg_start, reg_stop))
   bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = TRUE), what=c('rname', 'pos', 'qwidth', 'seq'), which=which)
   result <- Rsamtools::scanBam(bam_obj, param=bam_scan)[[1]]
   Rsamtools::close.BamFile(bam_obj)
   return(result)
}
