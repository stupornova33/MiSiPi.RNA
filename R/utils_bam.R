# close_bam is a wrapper function for Rsamtools close.BamFile function
# In rare instances, bam files haven't finished closing before they are trying to be deleted causing an error
# This is just an attempt to give the system a small amount of padding to finish closing a file before
# attempting to delete that file
# This doesn't feel like a very good way to deal with this error, but it will have to suffice for now
# Be aware that this could cause an infinite loop
# Consider adding a maximum number of iterations before exiting gracefully
.close_bam <- function(bam_obj) {
  if (Rsamtools::isOpen(bam_obj)) {
    Rsamtools::close.BamFile(bam_obj)
    while (Rsamtools::isOpen(bam_obj)) {
      Sys.sleep(0.5) # Sleep for half a second each time the file is still not closed when checking
    }
  }
  return()
}

# function to open bam file and verify that it was successfully accessed
# takes input file
# creates bam obj
# @param bamfile a BAM file
# @param logfile a string
# @return bam_obj

.open_bam <- function(bamfile, logfile){
  bam_obj <- Rsamtools::BamFile(bamfile)
  tryCatch(Rsamtools::open.BamFile(bam_obj), error = function(e) {
    msg <- conditionMessage(e)
    cat(msg, file = logfile, sep = "\n", append = TRUE)
    warning("Could not open Bamfile provided")
    stop(e)
  })
  return(bam_obj)
}
