# This function creates default null results to return in the event that no meaningful results could be created
#
# @return list

.null_hp_res <- function() {
  dicer_tbl <- data.frame(
    shift = seq(-4, 4),
    zscore = rep(0, times = 9),
    ml_zscore = rep(-33, times = 9)
  )
  
  hp_phased_tbl <- data.frame(
    phased_dist = seq(0, 50),
    phased_num = rep(0, times = 51),
    phased_z = rep(0, times = 51),
    phased_mlz = rep(-33, times = 51)
  )
  
  all_overlaps <- data.frame()
  
  null_results <- list(
    all_overlaps = all_overlaps,
    hp_overhangz = 0,
    hp_overhang_mlz = -33,
    hp_phasedz = 0,
    hp_phased_mlz = -33,
    hp_phased_tbl = hp_phased_tbl,
    phased_tbl.dist = hp_phased_tbl$phased_dist,
    phased_tbl.phased_z = hp_phased_tbl$phased_z,
    phased_tbl.phased_mlz = hp_phased_tbl$phased_mlz,
    dicer_tbl.shift = dicer_tbl$shift,
    dicer_tbl.zscore = dicer_tbl$zscore,
    dicer_tbl.ml_zscore = dicer_tbl$ml_zscore
  )
  
  return(null_results)
}
