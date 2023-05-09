#' function to open bam file and verify that it was successfully accessed
#' takes input file
#' creates bam obj
#' @param bamfile a BAM file
#' @param logfile
#' @return bam_obj



#' @export

OpenBamFile <- function(bamfile, logfile){
   bam_obj <- Rsamtools::BamFile(bamfile)
   Rsamtools::open.BamFile(bam_obj)
   if(Rsamtools::isOpen(bam_obj) == FALSE) {
      warning_msg <- paste("Warning", bamfile, "cannot be opened. Halting", sep = " ")
      cat(logfile, warning_msg, append = TRUE)
      return()
   } else {
      return(bam_obj)
   }
   Rsamtools::close.BamFile(bam_obj)
}

