# function to calculate the percent of reads that start with T
# takes forward and reverse data tables
#
# returns percent value
#
# @param forward_dt a data table of read data
# @param reverse_dt a data table of read data
# @return num

.first_nuc_T <- function(forward_dt, reverse_dt) {
  width <- NULL

  if (nrow(forward_dt) == 0) {
    all_dat <- reverse_dt
  } else if (nrow(reverse_dt) == 0) {
    all_dat <- forward_dt
  } else {
    all_dat <- rbind(forward_dt, reverse_dt) # %>% dplyr::filter(width >= 24 & width <= 30)
  }

  num_T <- all_dat %>% dplyr::filter(all_dat$first == "T")
  num <- nrow(num_T) / nrow(all_dat)
  return(num)
}
