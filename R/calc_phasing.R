#' calc phasing
#' plots output
#' takes two read dataframes
#' returns table of phasing scores
#'
#' @param df1 a dataframe of summarized reads with a duplicates column
#' @param df2 a dataframe of summarized reads with a duplicates column
#' @param n the number of bases to transform the read by
#' @return phased_counts

#' @export

calc_phasing <- function(df1, df2, n) {

  dist <- num.y <- num.x <- Zscore <- num <- NULL
  
  phased <- find_hp_overlaps(df1, df2, n)
  
  if (!nrow(phased) == 0) {
    print('making phased_counts')
    phased <- phased %>%
      dplyr::mutate(total_dupes = r1_dupes * r2_dupes) %>%
      dplyr::group_by(dist) %>%
      dplyr::summarize(num = sum(total_dupes)) %>%
      dplyr::mutate(dist = abs(dist))
  } else {
    print('setting phased_counts to zero')
    phased <- data.table::data.table(dist = 1, num = 0L)
  }
  
  #make the results data table
  print('making all_table')
  all_table <- data.table::data.table(dist = seq(0,50),
                                      num = rep(0, 51))
  
  if (nrow(phased) >= 1) {
    phased <- data.table::setDT(dplyr::full_join(phased,
                                                 all_table,
                                                 by = "dist", "num"))
    phased[is.na(phased)] <- 0
    phased <- phased %>%
      dplyr::select(-num.y)
    phased$Zscore <- calc_zscore(phased$num.x)
    phased <- phased %>%
      dplyr::rename(phased_dist = dist,
                    phased_num = num.x,
                    phased_z = Zscore)
  } else {
    phased <- data.table::data.table(dist=seq(0,50), num=rep(0, 51))
    phased$Zscore <- calc_zscore(phased$num)
    phased <- phased_counts %>%
      dplyr::rename(phased_dist = dist,
                    phased_num = num,
                    phased_z = Zscore)
  }
  
  return(phased)
}
