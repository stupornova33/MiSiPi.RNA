#' function to open bam file and verify that it was successfully accessed
#' takes input file
#' creates bam obj
#' @param bamfile a BAM file
#' @param logfile a string
#' @return bam_obj

#' @export

OpenBamFile <- function(bamfile, logfile){
   bam_obj <- Rsamtools::BamFile(bamfile)
   tryCatch(Rsamtools::open.BamFile(bam_obj), error = function(e) {
     msg <- conditionMessage(e)
     cat(msg, file = logfile, sep = "\n", append = TRUE)
     message("Could not open Bamfile provided")
     stop(e)
   })
   return(bam_obj)
}
