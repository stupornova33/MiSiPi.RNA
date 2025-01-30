# function to find hairpin overlaps
# takes two data tables of reads with start and stop positions
# converts to gr objects and finds overlaps, transforms positions back to original
# returns overlap table
#
# @param r1_dt a data frame containing chrom name, start, and stop
# @param r2_dt a data frame containing chrom name, start, and stop
# @param n the number of bases to transform the read by
# @return overlap table

.find_hp_overlaps <- function(r1_dt, r2_dt, n) {
  r1_stop <- r1_start <- dist <- start_r1 <- widthx <- start_r2 <- widthy <- NULL

  gr1_obj <- makeGRobj(r1_dt, "rname", "start", "end")
  gr2_obj <- makeGRobj(r2_dt, "rname", "start", "end")

  res <- GenomicRanges::findOverlaps(gr1_obj, gr2_obj,
    maxgap = -1L, minoverlap = 1L,
    type = c("any"),
    select = c("all"),
    ignore.strand = TRUE
  )

  # queryHits and subjectHits are functions to get the value at an index
  query_idx <- S4Vectors::queryHits(res)
  sub_idx <- S4Vectors::subjectHits(res)
  res <- NULL

  overlap_df <- data.frame(
    r1_start = IRanges::start(gr1_obj[query_idx]),
    r1_end = IRanges::end(gr1_obj[query_idx]) - n,
    r1_dupes = gr1_obj$n[query_idx],
    r2_start = IRanges::start(gr2_obj[sub_idx]),
    r2_width = IRanges::width(gr2_obj[sub_idx]),
    r2_end = IRanges::end(gr2_obj[sub_idx]),
    r2_dupes = gr2_obj$n[sub_idx]
  ) %>%
    dplyr::mutate(
      r1_width = r1_end - r1_start + 1,
      dist = r2_start - r1_end
    ) %>%
    dplyr::filter(
      dist >= 0 & dist < 50,
      r1_end <= r2_start
    )

  gr1_obj <- NULL
  gr2_obj <- NULL
  query_idx <- NULL
  sub_idx <- NULL

  return(overlap_df)
}
