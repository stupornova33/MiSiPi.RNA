#' function to calculate the read distribution
#' takes forward and reverse data frames
#' summarizes read size by count
#' returns table of size distribution
#' 
#' @param chrom_m a minus strand object
#' @param chrom_p a plus strand object
#' @return dizeDist



#' @export


get_read_dist <- function(chrom_m, chrom_p) {
   options(scipen = 999)
   
   plus_dt <- data.table::setDT(makeBamDF(chrom_p)) %>%
      base::subset(width <= 32 & width >= 18) %>% 
      dplyr::mutate(start = pos, end = pos + width - 1) %>%
      dplyr::select(-c(pos, seq))
   minus_dt <- data.table::setDT(makeBamDF(chrom_m)) %>%
      base::subset(width <= 32 & width >= 18) %>% 
      dplyr::mutate(start = pos, end = pos + width - 1) %>%
      dplyr::select(-c(pos, seq))
   
   all_dat <- rbind(plus_dt, minus_dt)
   sizeDist <- all_dat %>% dplyr::group_by(width, first) %>%
      dplyr::summarize(count = dplyr::n())
   
   widths <- c(rep(18, times=4), rep(19, times=4),rep(20, times=4),rep(21, times=4),rep(22, times=4),rep(23, times=4),rep(24, times=4),rep(25, times=4),
               rep(26, times=4),rep(27, times=4),rep(28, times=4),rep(29, times=4),rep(30, times=4),rep(31, times=4),rep(32, times=4))
   first <- c("A", "G", "C", "T")
   empty_dt <- data.frame(width = widths, first = rep(first, times = 15), count = 0)
   
   test <- merge(sizeDist, empty_dt, by = c('width', 'first'), all = TRUE) %>% dplyr::select(-c(count.y)) %>%
      dplyr::rename('count' = count.x)
   test["count"][is.na(test["count"])] <- 0
   
   size_dist <- test
   return(size_dist)
   
}
