#' null_mi_res
#'
#' creates a default object for return to miRNA_function in the event that no meaningful results could be calculated in a particular iteration
#'
#' @return default null result object
#' @export

null_mi_res <- function() {
  mfe <- 0
  z_df <- data.frame(overlap = seq(4,30), count = rep(0, times = 27))
  overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
  overhangs$z_score <- calc_zscore(overhangs$proper_count)
  return(list(list("mfe" = mfe, "perc_paired" = 0, "overhangs" = c(overhangs,z_df))))
}
