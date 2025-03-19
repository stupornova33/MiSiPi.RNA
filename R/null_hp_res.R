# This function creates default null results to return in the event that no meaningful results could be created
#
# @return list

.null_hp_res <- function() {
  null_results <- function() {
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
    
    return(list(
      MFE = 0,
      hp_overhangz = 0,
      hp_overhang_mlz = -33,
      hp_phasedz = 0,
      hp_phased_mlz = -33,
      phased_tbl.dist = phased_tbl$dist,
      phased_tbl.phased_z = phased_tbl$zscore,
      phased_tbl.phased_mlz = phased_tbl$ml_zscore,
      dicer_tbl.shift = dicer_tbl$shift,
      dicer_tbl.zscore = dicer_tbl$zscore,
      dicer_tbl.ml_zscore = dicer_tbl$ml_zscore,
      perc_paired = 0
    ))
  }

  minus_res <- null_results()
  plus_res <- null_results()
  
  return(list(
    minus_res = minus_res,
    plus_res = plus_res)
  )
}
