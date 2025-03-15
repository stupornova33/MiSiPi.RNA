# null_mi_res
#
# creates a default object for return to .miRNA function in the event that no meaningful results could be calculated in a particular iteration
#
# @return default null result object

.null_mi_res <- function() {
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
  
  return(list("mfe" = mfe, "perc_paired" = 0, "overhangs" = overhangs, "overlaps" = z_df))
}
