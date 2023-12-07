#' calc phasing
#' plots output
#' takes two read dataframes
#' returns table of phasing scores
#'
#' @param df1 a character passed in, either "+" or "-"
#' @param df2 a string
#' @return phased_counts

#' @export


calc_phasing <- function(df1, df2){
phased <- find_hp_overlaps(df1, df2)

dist <- start_r1 <- widthx <- start_r2 <- widthy <- num.y <- num.x <- Zscore <- num <- NULL

phased <- phased %>%
  dplyr::filter(dist >= 0) %>%
  dplyr::mutate(r1_end = r1_start + r1_width - 1, r2_end = r2_start + r2_width - 1)

if(!nrow(phased) == 0){
  print('making phased_counts')
  #filter_dt <- NULL
  phased_counts <- phased %>%
    dplyr::group_by(dist) %>%
    dplyr::summarize(num= dplyr::n()) %>%
    dplyr::mutate(dist = abs(dist)) %>%
    #dplyr::filter(dist <= 50)
    dplyr::filter(dist <= 63)


  #filter_r1_dt <- NULL
  #filter_r2_dt <- NULL
} else {
  print('setting phased_counts to zero')
  phased_counts <- data.table::data.table(dist = 1, num = 0L)
}


#make the results data table
print('making all_table')
all_table <- data.table::data.table(dist=seq(0,63), num=rep(0, 64))

if(nrow(phased_counts) >= 1){
  phased_counts <- data.table::setDT(dplyr::full_join(phased_counts, all_table, by = "dist", "num"))
  phased_counts[is.na(phased_counts)] <- 0
  phased_counts <- phased_counts %>% dplyr::select(-c(num.y))
  phased_counts$Zscore <- calc_zscore(phased_counts$num.x)
  phased_counts <- phased_counts %>% dplyr::rename(phased_dist = dist, phased_num = num.x, phased_z = Zscore)
} else {
  phased_counts <- data.table::data.table(dist=seq(0,63), num=rep(0, 64))
  phased_counts$Zscore <- calc_zscore(phased_counts$num)
  phased_counts <- phased_counts %>% dplyr::rename(phased_dist = dist, phased_num = num, phased_z = Zscore)
}

return(phased_counts)

}
