#' get_stranded_read_dist
#' takes forward and reverse data frames
#' summarizes read size by count
#' returns table of size distribution
#' @param bam_obj an object created by Rsamtools
#' @param chrom_name a string
#' @param reg_start an integer
#' @param reg_stop an integer
#' @return size_dist



#' @export


get_stranded_read_dist <- function(bam_obj, chrom_name, reg_start, reg_stop) {
  options(scipen = 999)
  
  chrom_p <- getChrPlus(bam_obj, chrom_name, reg_start, reg_stop)
  plus_dt <- data.table::setDT(makeBamDF(chrom_p)) %>%
    base::subset(width <= 32 & width >= 18) %>%
    dplyr::mutate(start = pos, end = pos + width - 1) %>%
    dplyr::select(-c(pos))
  chrom_p <- NULL
  
  chrom_m <- getChrMinus(bam_obj, chrom_name, reg_start, reg_stop)
  minus_dt <- data.table::setDT(makeBamDF(chrom_m)) %>%
    base::subset(width <= 32 & width >= 18) %>%
    dplyr::mutate(start = pos, end = pos + width - 1) %>%
    dplyr::select(-c(pos))
  chrom_m <- NULL
  
  
  plus_dist <- plus_dt %>% dplyr::group_by(width, first) %>% 
    dplyr::summarize(count = dplyr::n()) %>% dplyr::mutate(strand = "positive")
  
  minus_dist <- minus_dt %>% dplyr::group_by(width, first) %>% 
    dplyr::summarize(count = dplyr::n())
  
  minus_dist <- minus_dist %>% dplyr::mutate(neg_count = (count * -1), strand = "negative") %>%
    dplyr::select(-c(count)) %>% dplyr::rename("count" = "neg_count") 
  #plus_dt <- minus_dt <- NULL

  
 # all_dat <- rbind(plus_dist, minus_dist)

  widths <- c(rep(18, times=4), rep(19, times=4),rep(20, times=4),rep(21, times=4),rep(22, times=4),rep(23, times=4),rep(24, times=4),rep(25, times=4),
              rep(26, times=4),rep(27, times=4),rep(28, times=4),rep(29, times=4),rep(30, times=4),rep(31, times=4),rep(32, times=4))
  first <- c("A", "G", "C", "T")
  empty_dt <- data.frame(width = widths, first = rep(first, times = 15), count = 0)
  
  p_size_dist <- merge(plus_dist, empty_dt, by = c('width', 'first'), all = TRUE) %>%
    dplyr::select(-c(count.y)) %>%
    dplyr::rename('count' = count.x)
  p_size_dist["count"][is.na(p_size_dist["count"])] <- 0
  p_size_dist["strand"][is.na(p_size_dist["strand"])] <- "positive"
  
  m_size_dist <- merge(minus_dist, empty_dt, by = c('width', 'first'), all = TRUE) %>%
    dplyr::select(-c(count.y)) %>%
    dplyr::rename('count' = count.x)
  m_size_dist["count"][is.na(m_size_dist["count"])] <- 0
  m_size_dist["strand"][is.na(m_size_dist["strand"])] <- "negative"
  
  all_sizes <- rbind(p_size_dist, m_size_dist)
  
  
  
  sizeDist <- NULL
  
  return(all_sizes)
}
