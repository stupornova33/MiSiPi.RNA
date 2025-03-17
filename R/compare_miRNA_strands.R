.compare_miRNA_strands <- function(calling_func = c("miRNA", "all")) {
  
  calling_func <- match.arg(calling_func)
  output_dir <- switch(
    calling_func,
    "miRNA" = "miRNA_outputs",
    "all" = file.path("run_all", "miRNA_outputs")
  )
  
  plusFile <- file.path(output_dir, "miRNA_plus_dicerz.txt")
  minusFile <- file.path(output_dir, "miRNA_minus_dicerz.txt")
  
  pStrand <- FALSE
  mStrand <- FALSE
  
  if (file.exists(plusFile)) {
    plus_df <- read.delim(plusFile)
    pStrand = TRUE
  }
  
  if (file.exists(minusFile)) {
    minus_df <- read.delim(minusFile)
    mStrand = TRUE
  }
  
  if (!pStrand & !mStrand) {
    return()
  }
  
  output_file <- file.path(output_dir, "miRNA_dicerz.txt")
  
  
  if (!mStrand) {
    plus_df <- plus_df %>%
      dplyr::select(-count_avg)
    .write.quiet(plus_df, output_file)
    return()
  }
  
  if (!pStrand) {
    minus_df <- minus_df %>%
      dplyr::select(-count_avg)
    .write.quiet(minus_df, output_file)
    return()
  }
  
  original_loci <- unique(c(plus_df$original_locus, minus_df$original_locus))
  nRows <- length(original_loci)
  
  final_df <- data.frame(
    original_locus = character(nRows),
    most_abundant_locus = character(nRows),
    strand = character(nRows),
    `-4` = numeric(nRows),
    `-3` = numeric(nRows),
    `-2` = numeric(nRows),
    `-1` = numeric(nRows),
    `0` = numeric(nRows),
    `1` = numeric(nRows),
    `2` = numeric(nRows),
    `3` = numeric(nRows),
    `4` = numeric(nRows),
    check.names = FALSE
  )
  
  for (i in seq_along(original_loci)) {
    plus_row <- plus_df[plus_df$original_locus == original_loci[i], ]
    minus_row <- minus_df[minus_df$original_locus == original_loci[i], ]

    if (nrow(plus_row) == 0) {
      minus_row <- minus_row %>%
        dplyr::select(-count_avg)
      final_df[i,] <- minus_row
    } else if (nrow(minus_row) == 0) {
      plus_row <- plus_row %>%
        dplyr::select(-count_avg)
      final_df[i,] <- plus_row
    } else if (plus_row$count_avg > minus_row$count_avg) {
      plus_row <- plus_row %>%
        dplyr::select(-count_avg)
      final_df[i,] <- plus_row
    } else {
      minus_row <- minus_row %>%
        dplyr::select(-count_avg)
      final_df[i,] <- minus_row
    }
  }
  
  .write.quiet(final_df, output_file)
}
