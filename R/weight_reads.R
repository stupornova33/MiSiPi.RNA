# Takes a "collapsed" data.frame of reads summarized by count.
# Replicates them according to a factor provided by the user.
#
# @param filter_dt - Data Frame of summarized reads
# @param method - Character representing the method by which to expand the Data Frame
# @param norm - Integer chosen by the user which will be used to normalize the read counts.
# @param locus_length - Integer corresponding to the length of the region of interest
# @param locus_read_count - Integer
# @return expanded_dt A dataframe of reads which have been "uncollapsed" and replicated by the normalization factor.

.weight_reads <- function(filter_dt, weight_reads, locus_length, read_count) {
  if (missing(filter_dt)) {
    stop("Parameter `filter_dt` must be present")
  }
  if (missing(weight_reads)) {
    stop("Parameter `weight_reads` must be present")
  }

  if (is.numeric(weight_reads)) {
    method <- "value"
    norm <- weight_reads
  } else {
    method <- tolower(weight_reads)
  }

  if (method == "locus_norm") {
    if (missing(read_count)) {
      stop("Parameter `locus_read_count` must be present when using method `locus_norm`.")
    }
    stopifnot("Parameter `read_count` must be numeric." = is.numeric(read_count))
  }

  if (method == "value") {
    if (missing(locus_length)) {
      stop("Parameter `locus_length` must be present when weight_reads is user defined.")
    }
    stopifnot("locus_length must be numeric" = is.numeric(locus_length))
  }

  # output_msg <- switch(method,
  #   "none" = "No weighting of reads applied.",
  #   "locus_norm" = "Normalizing read count to locus.",
  #   "weight_by_prop" = "Weighting reads by proportion.",
  #   "value" = "User supplied custom weighting value for reads."
  # )
  # print(output_msg)

  # Weight reads
  expanded_dt <- filter_dt

  if (method == "value") {
    expanded_dt$rep_count <- round((expanded_dt$count * 10^9) / (norm * locus_length))
  } else if (method == "locus_norm") {
    options(scipen = 999)

    # CPM-like
    expanded_dt <- expanded_dt %>%
      dplyr::mutate(weighted_counts = round(count * (1 / locus_read_count) * 10^6)) %>%
      dplyr::select(-count) %>%
      dplyr::rename(count = weighted_counts)
  } else if (method == "weight_by_prop") {
    mn <- min(expanded_dt$count)
    mx <- max(expanded_dt$count)

    expanded_dt <- expanded_dt %>%
      dplyr::mutate(prop = count / mx) %>%
      dplyr::mutate(weighted_counts = prop * count) %>%
      dplyr::select(-count) %>%
      dplyr::rename(count = weighted_counts)

    # remove values < 1?
    expanded_dt <- expanded_dt[expanded_dt$count > 1, ]
  }

  expanded_dt <- expanded_dt %>%
    dplyr::arrange(count)
  
  expanded_dt <- rep_seq_reads(
    expanded_dt$count, expanded_dt$rname,
    expanded_dt$start, expanded_dt$end,
    expanded_dt$first, expanded_dt$seq
  ) %>%
    dplyr::mutate(width = end - start + 1)

  return(expanded_dt)
}
