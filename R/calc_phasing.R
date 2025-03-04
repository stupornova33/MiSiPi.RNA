# df1 - Data Frame of summarized reads with a duplicates column
# df2 - Data Frame of summarized reads with a duplicates column
# n - Integer representing number of bases to transform the read by
#
# returns Data Frame of phased scores

.calc_phasing <- function(df1, df2, n) {
  dist <- num.y <- num.x <- num <- NULL

  phased <- .find_hp_overlaps(df1, df2, n)

  # Convert the dupes columns to numeric so they can handle larger numbers
  # A few loci were giving combined duplicate counts into the trillions
  # The largest number R can store in an integer is system dependent,
  # but it is usually around 2 billion

  phased$r1_dupes <- as.numeric(phased$r1_dupes)
  phased$r2_dupes <- as.numeric(phased$r2_dupes)

  if (!nrow(phased) == 0) {
    phased <- phased %>%
      dplyr::mutate(total_dupes = r1_dupes * r2_dupes) %>%
      dplyr::group_by(dist) %>%
      dplyr::summarize(num = sum(total_dupes)) %>%
      dplyr::mutate(dist = abs(dist))
  } else {
    phased <- data.table::data.table(dist = 1, num = 0L)
  }

  # make the results data table
  all_table <- data.table::data.table(
    dist = seq(0, 50),
    num = rep(0, 51)
  )

  if (nrow(phased) >= 1) {
    phased <- data.table::setDT(dplyr::full_join(phased,
      all_table,
      by = "dist", "num"
    ))
    phased[is.na(phased)] <- 0
    phased <- phased %>%
      dplyr::select(-num.y)
    phased$phased_z <- .calc_zscore(phased$num.x)
    phased$phased_ml_z <- .calc_ml_zscore(phased$num.x)
    phased <- phased %>%
      dplyr::rename(
        phased_dist = dist,
        phased_num = num.x
      )
  } else {
    phased <- data.table::data.table(dist = seq(0, 50), num = rep(0, 51))
    phased$phased_z <- .calc_zscore(phased$num)
    phased$phased_ml_z <- .calc_ml_zscore(phased$num)
    phased <- phased_counts %>%
      dplyr::rename(
        phased_dist = dist,
        phased_num = num,
      )
  }

  return(phased)
}
