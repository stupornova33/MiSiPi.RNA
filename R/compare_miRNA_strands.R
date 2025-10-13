.compare_miRNA_strands <- function(output_dir) {
  
  plusFile <- file.path(output_dir, "miRNA", "miRNA_plus_dicerz.txt")
  minusFile <- file.path(output_dir, "miRNA", "miRNA_minus_dicerz.txt")
  
  plusFileExists <- FALSE
  minusFileExists <- FALSE
  
  if (file.exists(plusFile)) {
    plus_df <- read.delim(plusFile)
    plusFileExists = TRUE
  }
  
  if (file.exists(minusFile)) {
    minus_df <- read.delim(minusFile)
    minusFileExists = TRUE
  }
  
  if (!plusFileExists & !minusFileExists) {
    return()
  }
  
  output_file <- file.path(output_dir, "miRNA", "miRNA_dicerz.txt")
  
  if (!minusFileExists) {
    output_df <- plus_df %>%
      dplyr::select(-count_avg)
  } else if (!plusFileExists) {
    output_df <- minus_df %>%
      dplyr::select(-count_avg)
  } else if (plus_df$count_avg > minus_df$count_avg) {
    output_df <- plus_df %>%
      dplyr::select(-count_avg)
  } else {
    output_df <- minus_df %>%
      dplyr::select(-count_avg)
  }
  
  colnames(output_df) <- c("original_locus", "most_abundant_locus", "strand", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4")
  
  .write.quiet(output_df, output_file)
  
  # Remove this iteration's files
  if (plusFileExists) {
    file.remove(plusFile)
  }
  
  if (minusFileExists) {
    file.remove(minusFile)
  }
}
