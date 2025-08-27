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

.get_bam_scan <- function(bam_obj, chrom_name, reg_start, reg_stop, strand = c("minus", "plus", "-", "+")) {
  strand <- match.arg(strand)
  is_minus <- switch(strand,
                     "minus" = TRUE,
                     "plus" = FALSE,
                     "-" = TRUE,
                     "+" = FALSE
  )
  
  which <- GenomicRanges::GRanges(seqnames = chrom_name, IRanges::IRanges(reg_start, reg_stop))
  bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = is_minus, isUnmappedQuery = FALSE), what = c("rname", "pos", "qwidth", "seq"), which = which)
  result <- Rsamtools::scanBam(bam_obj, param = bam_scan)[[1]]
  return(result)
}

.get_filtered_bam_df <- function(bam_obj, chrom_name, reg_start, reg_stop, strand = c("minus", "plus", "-", "+"), min_width, max_width, include_seq = c(FALSE, TRUE)) {
  strand <- match.arg(strand)
  include_seq <- match.arg(include_seq)
  
  bam_scan <- .get_bam_scan(bam_obj, chrom_name, reg_start, reg_stop, strand)
  
  if (include_seq) {
    bam_df <- .get_bam_df_with_seq(bam_scan)
  } else {
    bam_df <- .get_bam_df(bam_scan)
  }
  
  bam_scan <- NULL
  
  bam_df <- bam_df %>%
    dplyr::filter(width >= min_width & width <= max_width) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1)
  
  return(bam_df)
}

# function to convert chrom obj into a data frame
# extracts rname, pos, width, and seq for siRNA dataframe
# Makes a dataframe out of the bam scan data
# @param chrom_obj an object created from a bam obj
# @return bam_df

.get_bam_df_with_seq <- function(chrom_obj) {
  bam_df <- .get_bam_df(chrom_obj)
  bam_df <- bam_df %>%
    dplyr::mutate(seq = chrom_obj$seq)
  
  return(bam_df)
}

# function to convert chrom obj into a data frame
# extracts rname, pos, width, seq, and the first letter of the sequence
# Makes a dataframe out of the bam scan data
# @param chrom_obj an object created from a bam obj
# @return bam_df

.get_bam_df <- function(chrom_obj) {
  # get just first nuc of sequence to save memory
  s1 <- substr(chrom_obj$seq, 1, 1)
  bam_df <- data.frame(
    "rname" = as.character(chrom_obj$rname),
    "pos" = chrom_obj$pos,
    "width" = chrom_obj$qwidth,
    "first" = s1
  )
  return(bam_df)
}
