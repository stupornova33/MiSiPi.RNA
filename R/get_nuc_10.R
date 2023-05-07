#' function to get percentage of reads with A at 10th nucleotide position
#' takes reverse dt and forward dt
#' returns percentage
#' 
#' @param forward_dt a data table of read data
#' @param reverse_dt a data table of read data
#' @return percentage



#' @export

get_nuc_10 <- function(forward_dt, reverse_dt) {
   all_dat <- rbind(forward_dt, reverse_dt) %>% dplyr::filter(width >= 24 & width <= 30)
   total <- 0
   accumulate_nucs <- function(datseq) {
      nuc <- stringr::str_split(datseq, pattern = "")[[1]][10]
      if (nuc == "A") total <<- total + 1
   }
   lapply(all_dat$seq, accumulate_nucs)
   return(total/nrow(all_dat))
}
