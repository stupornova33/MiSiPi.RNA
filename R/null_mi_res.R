# null_mi_res
#
# creates a default object for return to .miRNA function in the event that no meaningful results could be calculated in a particular iteration
#
# @return default null result object

.null_mi_res <- function(region_str, strand, wkdir) {
  mfe <- 0
  z_df <- data.frame(Overlap = seq(4, 30), count = rep(0, times = 27))
  z_df$zscore <- .calc_zscore(z_df$count)
  z_df$ml_zscore <- .calc_ml_zscore(z_df$count)
  
  overhangs <- data.frame(
    shift = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
    proper_count = c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    improper_count = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
  )
  overhangs$zscore <- .calc_zscore(overhangs$proper_count)
  overhangs$ml_zscore <- .calc_ml_zscore(overhangs$proper_count)
  
  # Write null results to files
  #### Overlap Probability ####
  z_output <- z_df %>%
    dplyr::select(c(-ml_zscore, count))
  
  zdf_output <- as.data.frame(t(z_output[,-1]))
  colnames(zdf_output) <- z_output$Overlap
  
  zdf_file <- file.path(wkdir, "overlap_probability.txt")
  .write.quiet(zdf_output, zdf_file)
  
  #### Dicer Zscore ####
  overhang_output <- data.frame(t(overhangs$ml_zscore))
  colnames(overhang_output) <- overhangs$shift
  overhang_output$original_locus <- region_str
  overhang_output$most_abundant_locus <- NA
  overhang_output$strand <- strand
  overhang_output$count_avg <- 0
  overhang_output <- overhang_output[, c(10, 11, 12, 13, 1:9)]
  
  dice_file <- switch(
    strand,
    "+" = "miRNA_plus_dicerz.txt",
    "-" = "miRNA_minus_dicerz.txt"
  )
  dice_file <- file.path(wkdir, dice_file)
  .write.quiet(overhang_output, dice_file)
  
  return(list("mfe" = mfe, "perc_paired" = 0, "overhangs" = overhangs, "overlaps" = z_df))
}
