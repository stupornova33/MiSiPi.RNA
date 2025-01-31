# function to convert chrom obj into a data frame
# extracts rname, pos, width, and seq for siRNA dataframe
# Makes a dataframe out of the bam scan data
# @param chrom_obj an object created from a bam obj
# @return bam_df

.make_si_BamDF <- function(chrom_obj) {
  # get just first nuc of sequence to save memory
  s1 <- substr(chrom_obj$seq, 1, 1)
  bam_df <- data.frame("rname" = chrom_obj$rname, "pos" = chrom_obj$pos, "width" = chrom_obj$qwidth, "seq" = chrom_obj$seq, "first" = s1)
  return(bam_df)
}
