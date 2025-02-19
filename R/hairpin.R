# run the hairpin function
# runs hairpin function on both strands
# takes a gtf file, bam file, and genome file specified by user
#
# @param vars a list

hairpin <- function(vars) {
  wkdir <- "hairpin_outputs/"
  if (!dir.exists(wkdir) == TRUE) dir.create(wkdir)

  logfile <- "hairpin_logfile.txt"
  if (!file.exists(logfile) == TRUE) file.create(paste0(wkdir, logfile))

  invisible(
    mapply(
      dual_strand_hairpin, vars[[1]], vars[[2]], vars[[3]], vars[[4]], vars[[9]], vars[[10]], logfile, wkdir, vars[[6]],
      vars[[7]], vars[[14]], vars[[15]], vars[[18]]
    )
  )
}
