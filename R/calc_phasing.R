#' calc phasing
#' plots output
#' takes two read dataframes
#' returns table of phasing scores
#'
#' @param df1 a character passed in, either "+" or "-"
#' @param df2 a string
#' @param n the number of bases to transform the read by
#' @return phased_counts

#' @export


calc_phasing <- function(df1, df2, n){


dist <- start_r1 <- widthx <- start_r2 <- widthy <- num.y <- num.x <- z_score <- num <- NULL

# check if number of reads is > 3000 and if so, subsample
df1 <- df1[sample(1:nrow(df1)),]
df2 <- df2[sample(1:nrow(df2)),]
df1 <- utils::head(df1, 4000)
df2 <- utils::head(df2, 4000)
phased <- find_hp_overlaps(df1, df2, n)

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
    dplyr::filter(dist <= 50)


  #filter_r1_dt <- NULL
  #filter_r2_dt <- NULL
} else {
  print('setting phased_counts to zero')
  phased_counts <- data.table::data.table(dist = 1, num = 0L)
}


#make the results data table
print('making all_table')
all_table <- data.table::data.table(dist=seq(0,50), num=rep(0, 51))

if(nrow(phased_counts) >= 1){
  phased_counts <- data.table::setDT(dplyr::full_join(phased_counts, all_table, by = "dist", "num"))
  phased_counts[is.na(phased_counts)] <- 0
  phased_counts <- phased_counts %>% dplyr::select(-c(num.y))
  phased_counts$z_score <- calc_zscore(phased_counts$num.x)
  phased_counts <- phased_counts %>% dplyr::rename(phased_dist = dist, phased_num = num.x, phased_z = z_score)
} else {
  phased_counts <- data.table::data.table(dist=seq(0,50), num=rep(0, 51))
  phased_counts$z_score <- calc_zscore(phased_counts$num)
  phased_counts <- phased_counts %>% dplyr::rename(phased_dist = dist, phased_num = num, phased_z = z_score)
}

return(phased_counts)

}
