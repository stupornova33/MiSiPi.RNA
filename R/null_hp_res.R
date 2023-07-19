#' @return list

#' @export

null_hp_res <- function(){
  neg_results <- function(){
    MFE <- 0
    minus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
    minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
    #minus_phased_tbl <- data.table::data.table(phased_dist = seq(1,50), phased_num = rep(0,50), phased_z = rep("NaN",50))
    #minus_hp_phased_counts <- sum(minus_phased_tbl$phased_num[1:4])
    minus_hp_phased_z <- -33
    minus_hp_overhangz <- -33
    return(c(minusMFE = MFE, minus_overhangz = minus_hp_overhangz, minus_hp_phased_z = minus_hp_phased_z, perc_paired = 0))
  }
  
  pos_results <- function(){
    MFE <- 0
    plus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
    plus_overhangs$zscore <- calc_zscore(plus_overhangs$proper_count)
    #plus_phased_tbl <- data.table::data.table(phased_dist = seq(1,50), phased_num = rep(0,50), phased_z = rep("NaN",50))
    #plus_hp_phased_counts <- sum(minus_phased_tbl$phased_num[1:4])
    plus_hp_phased_z <- -33
    plus_hp_overhangz <- -33
    return(c(plusMFE = MFE, plus_hp_overhangz = plus_hp_overhangz, plus_hp_phased_z = plus_hp_phased_z, perc_paired = 0))
  }
  
  neg_res <- neg_results()
  pos_res <- pos_results()
  return(list(neg_res, pos_res))
  
}