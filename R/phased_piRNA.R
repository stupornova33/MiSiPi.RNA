# phased_piRNA
# calls phased_piRNA function on both strand, creates output wkdir and logfile
# @param vars a list

phased_piRNA <- function(vars) {
  logfile <- "phased_piRNA_logfile.txt"

  wkdir <- "phased_piRNA_outputs/"
  if (!dir.exists(wkdir) == TRUE) dir.create(wkdir)

  logfile <- "phased_piRNA_logfile.txt"
  if (!file.exists(logfile) == TRUE) file.create(paste0(wkdir, logfile))

  invisible(mapply(.phased_piRNA, "+", vars[[1]], vars[[2]], vars[[3]], vars[[10]], logfile, wkdir, vars[[6]]))
  invisible(mapply(.phased_piRNA, "-", vars[[1]], vars[[2]], vars[[3]], vars[[10]], logfile, wkdir, vars[[6]]))
}
