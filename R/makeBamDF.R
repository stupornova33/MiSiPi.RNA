#' function to convert chrom obj into a data frame 
#' extracts rname, pos, width, seq, and the first letter of the sequence
#' finds max row, column, and count of highest cell
#' 
#' @param chrom_obj an object created from a bam obj
#' @return bam_df

#' @export

makeBamDF <- function(chrom_obj) { # Makes a dataframe out of the bam scan data
   # get just first nuc of sequence to save memory 
   s1 <- substr(chrom_obj$seq, 1,1)
   bam_df <- data.frame('rname' = chrom_obj$rname,'pos'= chrom_obj$pos, 'width' = chrom_obj$qwidth,'seq' = chrom_obj$seq, 'first' = s1)
   return(bam_df)
}
