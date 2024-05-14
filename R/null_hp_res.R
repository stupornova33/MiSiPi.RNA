#' This function creates default null results to return in the event that no meaningful results could be created
#'
#' @return list
#' @export

null_hp_res <- function(){
  neg_results <- function(){
    MFE <- 0
    dicer_tbl <- data.frame(shift = c(seq(-4,4)), zscore = c(rep(NA, times = 9)))
    phased_tbl <- data.frame(phased_dist = c(seq(0,50)), phased_num = c(rep(0, times = 51)), zscore = c(rep(NA, times = 51)))
    #phased_tbl.phased_z <- c(rep("NaN", times = 51))

    minus_hp_phasedz <- -33
    minus_hp_overhangz <- -33
    return(list(minusMFE = MFE, minus_hp_overhangz = minus_hp_overhangz, minus_hp_phasedz = minus_hp_phasedz,
             dicer_tbl.shift = dicer_tbl$shift, dicer_tbl.zscore = dicer_tbl$zscore, phased_tbl.phased_z = phased_tbl$zscore,
             phased_tbl.dist = phased_tbl$dist, perc_paired = 0))
  }



  pos_results <- function(){
    MFE <- 0
    dicer_tbl <- data.frame(shift = c(seq(-4,4)), zscore = c(rep(NA, times = 9)))
    phased_tbl <- data.frame(dist = c(seq(0,50)), phased_num = c(rep(0, times = 51)), zscore = c(rep(NA, times = 51)))
    plus_hp_phasedz <- -33
    plus_hp_overhangz <- -33
    return(list(plusMFE = MFE, plus_hp_overhangz = plus_hp_overhangz, plus_hp_phasedz = plus_hp_phasedz,
                dicer_tbl.shift = dicer_tbl$shift, dicer_tbl.zscore = dicer_tbl$zscore, phased_tbl.phased_z = phased_tbl$zscore,
                phased_tbl.dist = phased_tbl$dist, perc_paired = 0))
  }

  neg_res <- neg_results()
  pos_res <- pos_results()
  return(list(neg_res, pos_res))

}
