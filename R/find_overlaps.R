# function to find regular overlaps
# takes two data tables of reads with start and stop positions
# converts to gr objects and finds overlaps
# returns overlap table
#
# @param r1_dt a data frame which contains a read_name, start, and stop
# @param r2_dt a data frame which contains a read_name, start, and stop
# @return overlap table
# @export

.find_overlaps <- function(r1_dt, r2_dt = NULL) {
  # Currently it is expecting the r2_dt to have a column of summarized counts called "n"
  # So all duplicate rows are removed. The output of this function keeps that column in
  # So it can be used downstream
  # TODO We either need to incorporate a duplicate column for both dt's OR
  # we need to check ahead of time which dt has more duplicates, and only summarize that one
  # if we go the second route, then we'll need another variable to track which dt is summarized

  # Attempting to work with both DTs as summarized

  if (!is.null(r2_dt)) {
    gr1_obj <- .makeGRobj(r1_dt, "rname", "start", "end")
    gr2_obj <- .makeGRobj(r2_dt, "rname", "start", "end")

    res <- GenomicRanges::findOverlaps(gr1_obj, gr2_obj,
      maxgap = -1L, minoverlap = 4L,
      type = c("any"),
      select = c("all"),
      ignore.strand = TRUE
    )
  } else {
    gr1_obj <- .makeGRobj(r1_dt, "rname", "start", "end")
    res <- GenomicRanges::findOverlaps(gr1_obj,
      maxgap = -1L, minoverlap = 4L,
      type = c("any"), select = c("all"),
      ignore.strand = TRUE,
      drop.self = TRUE,
      drop.redundant = TRUE
    )
    gr2_obj <- gr1_obj
  }

  # queryHits and subjectHits are functions to get the value at an index
  query_idx <- S4Vectors::queryHits(res)
  sub_idx <- S4Vectors::subjectHits(res)
  res <- NULL

  overlap_df <- data.frame(
    r1_start = IRanges::start(gr1_obj[query_idx]),
    r1_width = IRanges::width(gr1_obj[query_idx]),
    r1_end = IRanges::end(gr1_obj[query_idx]),
    r1_dupes = gr1_obj$n[query_idx],
    r2_start = IRanges::start(gr2_obj[sub_idx]),
    r2_width = IRanges::width(gr2_obj[sub_idx]),
    r2_end = IRanges::end(gr2_obj[sub_idx]),
    r2_dupes = gr2_obj$n[sub_idx]
  )

  gr1_obj <- NULL
  gr2_obj <- NULL
  query_idx <- NULL
  sub_idx <- NULL

  return(overlap_df)
}
