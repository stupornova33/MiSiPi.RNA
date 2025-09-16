# Wrapper function to handle timeout errors of .read_densityBySize

.read_density_by_size <- function(chrom_name, reg_start, reg_stop, input_file, wkdir, logfile, timeout = 3600) {
  
  # This function will be called in the finally block to cleanup any remains of the timed out function
  cleanup_bams <- function() {
    bam_path <- file.path(wkdir, "tmp_bam")
    
    # Gather and close bam files in case they are still open
    bam_files <- list.files(path = bam_path, pattern = "\\.bam$")
    
    #for (bam_file in bam_files) {
    #  if (Rsamtools::isOpen(bam_file)) {
    #    .close_bam(bam_file)
    #  }
    #}
    
    # Add index files to list and delete files and directory
    bam_files <- append(bam_files, list.files(path = bam_path, pattern = "\\.bai$"))
    
    unlink(file.path(bam_path, bam_files), force = TRUE)
    fs::dir_delete(bam_path)
  }
  
  
  # Set a default result in case read_densityBySize times out and fails to return anything
  density_data <- NULL
  
  tryCatch(
    density_data <- R.utils::withTimeout(
      .read_densityBySize(chrom_name, reg_start, reg_stop, input_file, wkdir, logfile),
      timeout = timeout
    ),
    # Log the errors and move on
    error = function(e) {
      error_msg <- paste("ERROR: Timeout of", timeout, "seconds was exceeded by .read_density_by_size. Timeout can be increased in set_vars.\nSetting null plot for read density and continuing.")
      cat(file = logfile, error_msg)
    },
    # Some temporary bam files might still exist that need to be removed in the finally block
    finally = cleanup_bams()
  )
  
  # Will either be NULL or a data frame
  # Plot density_data
  if (is.null(density_data)) {
    density_plot <- null_plot("Read Density Plot", "Timeout exceeded during processing")
  } else {
    if ((reg_stop - reg_start) > 7000) {
      density_plot <- .plot_large_density(density_data, reg_start, reg_stop)
    } else {
      density_plot <- .plot_density(density_data, reg_start, reg_stop)
    }
  }
  
  density_data <- NULL
  
  return(density_plot)
}

# function to filter reads by size and plot pileup density
# @param chrom_name a string
# @param reg_start a string
# @param reg_stop a string
# @param input_file a string
# @param wkdir a string
# @param logfile a string
# @return all_df

.read_densityBySize <- function(chrom_name, reg_start, reg_stop, input_file, wkdir, logfile) {
  pos <- width <- rname <- NULL
  
  filter_bamfile <- function(input_file, size1, size2, strand) {
    seqnames <- NULL
    which <- GenomicRanges::GRanges(
      seqnames = chrom_name,
      IRanges::IRanges(reg_start, reg_stop)
    )
    
    # Consider finding an alternative to FilterRules
    # FilterRules loads all reads into memory before filtering
    filters <- S4Vectors::FilterRules(list(MinWidth = function(x) {
      (x$qwidth >= size1 & x$qwidth <= size2)
    }))

    if (strand == "+") {
      filename <- paste(size1, size2, "pos.bam", sep = "_")
      filepath <- file.path(bam_path, filename)

      Rsamtools::filterBam(
        file = input_file,
        destination = filepath,
        filter = filters,
        param = Rsamtools::ScanBamParam(
          flag = Rsamtools::scanBamFlag(isMinusStrand = FALSE),
          what = c("pos", "qwidth"),
          which = which
        )
      )
    } else {
      filename <- paste(size1, size2, "neg.bam", sep = "_")
      filepath <- file.path(bam_path, filename)

      Rsamtools::filterBam(
        file = input_file,
        destination = filepath,
        filter = filters,
        param = Rsamtools::ScanBamParam(
          flag = Rsamtools::scanBamFlag(isMinusStrand = TRUE),
          what = c("pos", "qwidth"),
          which = which
        )
      )
    }

    new_bam_obj <- .open_bam(filepath, logfile)
    return(new_bam_obj)
  }

  # Create a subdirectory for temp bam files
  bam_path <- file.path(wkdir, "tmp_bam")
  if (!dir.exists(bam_path)) dir.create(bam_path)

  filtered_pos_18_19_bam <- filter_bamfile(input_file, 18, 19, "+")
  filtered_pos_20_22_bam <- filter_bamfile(input_file, 20, 22, "+")
  filtered_pos_23_25_bam <- filter_bamfile(input_file, 23, 25, "+")
  filtered_pos_26_32_bam <- filter_bamfile(input_file, 26, 32, "+")

  filtered_neg_18_19_bam <- filter_bamfile(input_file, 18, 19, "-")
  filtered_neg_20_22_bam <- filter_bamfile(input_file, 20, 22, "-")
  filtered_neg_23_25_bam <- filter_bamfile(input_file, 23, 25, "-")
  filtered_neg_26_32_bam <- filter_bamfile(input_file, 26, 32, "-")


  make_bam_pileup <- function(bam, strand) {
    seqnames <- pos <- count <- NULL

    which <- GenomicRanges::GRanges(seqnames = chrom_name, IRanges::IRanges(reg_start, reg_stop))
    bam_scan <- Rsamtools::ScanBamParam(
      flag = Rsamtools::scanBamFlag(isMinusStrand = strand == "-"),
      what = c("pos"),
      which = which
    )

    params <- Rsamtools::PileupParam(
      max_depth = 4000,
      min_base_quality = 20,
      min_mapq = 0,
      min_nucleotide_depth = 1,
      distinguish_strands = TRUE,
      distinguish_nucleotides = TRUE,
      ignore_query_Ns = TRUE,
      include_deletions = TRUE,
      include_insertions = FALSE,
      left_bins = NULL,
      query_bins = NULL,
      cycle_bins = NULL
    )

    pileups <- Rsamtools::pileup(bam,
      index = (stringr::str_c(input_file, "", ".bai")),
      scanBamParam = bam_scan,
      pileupParam = params
    ) %>%
      dplyr::select(-c(seqnames, strand))

    # Summarize duplicate positions including minor alleles
    dt <- pileups %>%
      dplyr::group_by(pos) %>%
      dplyr::summarise(count = sum(count))
    return(dt)
  }

  ## 18-19 nt
  pos_18_19_pileup <- make_bam_pileup(filtered_pos_18_19_bam, "+")
  neg_18_19_pileup <- make_bam_pileup(filtered_neg_18_19_bam, "-")

  ## 20-22 nt
  pos_20_22_pileup <- make_bam_pileup(filtered_pos_20_22_bam, "+")
  neg_20_22_pileup <- make_bam_pileup(filtered_neg_20_22_bam, "-")

  ## 23-25nt
  pos_23_25_pileup <- make_bam_pileup(filtered_pos_23_25_bam, "+")
  neg_23_25_pileup <- make_bam_pileup(filtered_neg_23_25_bam, "_")

  ## 26-32
  pos_26_32_pileup <- make_bam_pileup(filtered_pos_26_32_bam, "+")
  neg_26_32_pileup <- make_bam_pileup(filtered_neg_26_32_bam, "-")


  # close the bamfiles to be deleted
  .close_bam(filtered_pos_18_19_bam)
  .close_bam(filtered_neg_18_19_bam)
  .close_bam(filtered_pos_20_22_bam)
  .close_bam(filtered_neg_20_22_bam)
  .close_bam(filtered_pos_23_25_bam)
  .close_bam(filtered_neg_23_25_bam)
  .close_bam(filtered_pos_26_32_bam)
  .close_bam(filtered_neg_26_32_bam)
  
  # Clean up, clean up
  del_files <- list.files(path = bam_path, pattern = "\\.bam$")
  del_files <- append(del_files, list.files(path = bam_path, pattern = "\\.bai$"))

  unlink(file.path(bam_path, del_files), force = TRUE)
  fs::dir_delete(bam_path)

  empty_dat <- data.frame(pos = c(seq(reg_start, reg_stop)))

  merge_neg_table <- function(empty_dat, res_dat) {
    count <- NULL
    neg_res <- merge(empty_dat, res_dat, by = "pos", all.x = TRUE) %>%
      dplyr::mutate(count = count * -1) %>%
      dplyr::mutate(size = "all")

    neg_res["count"][is.na(neg_res["count"])] <- 0
    return(neg_res)
  }

  merge_pos_table <- function(empty_dat, res_dat) {
    pos_res <- merge(empty_dat, res_dat, by = "pos", all.x = TRUE) %>%
      dplyr::mutate(size = "all")
    pos_res["count"][is.na(pos_res["count"])] <- 0

    return(pos_res)
  }

  pos_18_19_res <- merge_pos_table(empty_dat, pos_18_19_pileup)
  pos_18_19_res$size <- "18_19"
  neg_18_19_res <- merge_neg_table(empty_dat, neg_18_19_pileup)
  neg_18_19_res$size <- "18_19"
  pos_20_22_res <- merge_pos_table(empty_dat, pos_20_22_pileup)
  pos_20_22_res$size <- "20_22"
  neg_20_22_res <- merge_neg_table(empty_dat, neg_20_22_pileup)
  neg_20_22_res$size <- "20_22"
  pos_23_25_res <- merge_pos_table(empty_dat, pos_23_25_pileup)
  pos_23_25_res$size <- "23_25"
  neg_23_25_res <- merge_neg_table(empty_dat, neg_23_25_pileup)
  neg_23_25_res$size <- "23_25"
  pos_26_32_res <- merge_pos_table(empty_dat, pos_26_32_pileup)
  pos_26_32_res$size <- "26_32"
  neg_26_32_res <- merge_neg_table(empty_dat, neg_26_32_pileup)
  neg_26_32_res$size <- "26_32"

  pos_18_19_pileup <- neg_18_19_pileup <- pos_20_22_pileup <- neg_20_22_pileup <- NULL
  pos_23_25_pileup <- neg_23_25_pileup <- pos_26_32_pileup <- neg_26_32_pileup <- NULL

  pos_counts <- dplyr::bind_rows(pos_18_19_res, pos_20_22_res, pos_23_25_res, pos_26_32_res)
  neg_counts <- dplyr::bind_rows(neg_18_19_res, neg_20_22_res, neg_23_25_res, neg_26_32_res)

  pos_18_19_res <- neg_18_19_res <- pos_20_22_res <- neg_20_22_res <- NULL
  pos_23_25_res <- neg_23_25_res <- pos_26_32_res <- neg_26_32_res <- NULL

  pos_counts$strand <- "pos"
  neg_counts$strand <- "neg"

  df <- dplyr::bind_rows(pos_counts, neg_counts)

  return(df)
}
