#' function to calculate the percent of reads that start with T
#' takes forward and reverse data tables
#' 
#' returns percent value
#' 
#' @param forward_dt a data table of read data
#' @param reverse_dt a data table of read data
#' @return num



#' @export


first_nuc_T <- function(forward_dt, reverse_dt){
   width <- NULL
   all_dat <- rbind(forward_dt, reverse_dt) %>% dplyr::filter(width >= 24 & width <= 30)
   num_T <- all_dat %>% dplyr::filter(all_dat$first == "T")
   num <- nrow(num_T)/nrow(all_dat)
   return(num)
}
