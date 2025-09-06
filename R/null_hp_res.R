# This function creates default null results to return in the event that no meaningful results could be created
#
# @return list

.null_hp_res <- function() {
  dicer_tbl <- data.frame(
    shift = seq(-4, 4),
    zscore = rep(0, times = 9),
    ml_zscore = rep(-33, times = 9)
  )
  
  phased_tbl <- data.frame(
    dist = seq(0, 50),
    phased_num = rep(0, times = 51),
    zscore = rep(0, times = 51),
    ml_zscore = rep(-33, times = 51)
  )
  
  all_overlaps <- data.frame()
  
  null_results <- list(
    all_overlaps = all_overlaps,
    hp_overhangz = 0,
    hp_overhang_mlz = -33,
    hp_phasedz = 0,
    hp_phased_mlz = -33,
    phased_tbl.dist = phased_tbl$dist,
    phased_tbl.phased_z = phased_tbl$zscore,
    phased_tbl.phased_mlz = phased_tbl$ml_zscore,
    dicer_tbl.shift = dicer_tbl$shift,
    dicer_tbl.zscore = dicer_tbl$zscore,
    dicer_tbl.ml_zscore = dicer_tbl$ml_zscore
  )
  
  return(null_results)
}
