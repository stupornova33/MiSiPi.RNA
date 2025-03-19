.compare_miRNA_strands <- function(chrom, start, stop, calling_func = c("miRNA", "all")) {
  
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
  
  # Create a vector of region strings for iterating through
  original_loci <- unique(c(plus_df$original_locus, minus_df$original_locus))
  bed_loci <- .get_region_string(chrom, start, stop)
  nRows <- length(bed_loci)
  
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
  
  for (locus in original_loci) {
    plus_row <- plus_df[plus_df$original_locus == locus, ]
    minus_row <- minus_df[minus_df$original_locus == locus, ]

    nRow_plus <- nrow(plus_row)
    nRow_minus <- nrow(minus_row)
    
    # Get each row index where this locus is present in the bed file
    bed_idx <- which(bed_loci %in% locus)
    
    # If we have more than 1 result, then they are duplicates
    # Take the first row
    plus_row <- plus_row %>%
      dplyr::slice_head(n = 1)
    minus_row <- minus_row %>%
      dplyr::slice_head(n = 1)
    
    # There will never be a situation where both nRows are 0
    if (nRow_plus == 0) {
      minus_row <- minus_row %>%
        dplyr::select(-count_avg)

      final_df[bed_idx,] <- minus_row
    } else if (nRow_minus == 0) {
      plus_row <- plus_row %>%
        dplyr::select(-count_avg)
      
      final_df[bed_idx,] <- plus_row
    } else if (plus_row$count_avg > minus_row$count_avg) {
      plus_row <- plus_row %>%
        dplyr::select(-count_avg)
      final_df[bed_idx,] <- plus_row
    } else {
      minus_row <- minus_row %>%
        dplyr::select(-count_avg)
      final_df[bed_idx,] <- minus_row
    }
  }
  
  .write.quiet(final_df, output_file)
}
