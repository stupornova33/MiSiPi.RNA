#' get_read_size_dist
#' takes forward and reverse data frames
#' summarizes read size by count
#' returns table of size distribution
#' @param f_df a data.frame of forward reads
#' @param r_df a data.frame of reverse reads
#' @return size_dist

#' @export

get_read_size_dist <- function(f_df, r_df) {
  options(scipen = 999)
  count.y <- count.x <- NULL
  
  sizeDist <- dplyr::bind_rows(f_df, r_df) %>%
    dplyr::group_by(width, first) %>%
    dplyr::summarize(count = dplyr::n())
  
  widths <- c(rep(18, times=4), rep(19, times=4),rep(20, times=4),rep(21, times=4),rep(22, times=4),rep(23, times=4),rep(24, times=4),rep(25, times=4),
              rep(26, times=4),rep(27, times=4),rep(28, times=4),rep(29, times=4),rep(30, times=4),rep(31, times=4),rep(32, times=4))
  first <- c("A", "G", "C", "T")
  empty_dt <- data.frame(width = widths, first = rep(first, times = 15), count = 0)
  
  size_dist <- merge(sizeDist, empty_dt, by = c('width', 'first'), all = TRUE) %>%
    dplyr::select(-c(count.y)) %>%
    dplyr::rename('count' = count.x)
  size_dist["count"][is.na(size_dist["count"])] <- 0
  
  sizeDist <- NULL
  
  return(size_dist)
}
