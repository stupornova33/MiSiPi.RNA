# get_read_dist
# takes forward and reverse data frames
# summarizes read size by count
# returns table of size distribution
# @param bam_obj an object created by Rsamtools
# @param chrom_name a string
# @param reg_start an integer
# @param reg_stop an integer
# @return size_dist

.get_read_dist <- function(bam_obj, chrom_name, reg_start, reg_stop) {
  options(scipen = 999)
  width <- pos <- count.y <- count.x <- NULL

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

  all_dat <- rbind(plus_dt, minus_dt)

  plus_dt <- minus_dt <- NULL

  sizeDist <- all_dat %>%
    dplyr::group_by(width, first) %>%
    dplyr::summarize(count = dplyr::n())

  all_dat <- NULL

  widths <- c(
    rep(18, times = 4), rep(19, times = 4), rep(20, times = 4), rep(21, times = 4), rep(22, times = 4), rep(23, times = 4), rep(24, times = 4), rep(25, times = 4),
    rep(26, times = 4), rep(27, times = 4), rep(28, times = 4), rep(29, times = 4), rep(30, times = 4), rep(31, times = 4), rep(32, times = 4)
  )
  first <- c("A", "G", "C", "T")
  empty_dt <- data.frame(width = widths, first = rep(first, times = 15), count = 0)

  size_dist <- merge(sizeDist, empty_dt, by = c("width", "first"), all = TRUE) %>%
    dplyr::select(-c(count.y)) %>%
    dplyr::rename("count" = count.x)
  size_dist["count"][is.na(size_dist["count"])] <- 0

  sizeDist <- NULL

  return(size_dist)
}
