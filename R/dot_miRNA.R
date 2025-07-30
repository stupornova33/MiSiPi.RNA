# .miRNA
# processes reads according to miRNA_algorithm, new with rnaplot function
# returns plots
# @param chrom_name a string
# @param reg_start a whole number
# @param reg_stop a whole number
# @param length length of chromosome of interest
# @param strand a character passed in, "+" or "-"
# @param genome_file a fasta file of chrom sequences
# @param bam_file a BAM file
# @param logfile a string
# @param wkdir a string
# @param plot_output a bool, default = TRUE
# @param path_to_RNAfold a string
# @param path_to_RNAplot a string
# @param write_fastas a bool, TRUE or FALSE
# @param weight_reads Determines whether read counts will be weighted and with which method. Valid options are "weight_by_prop", "locus_norm", a user-defined value, or "none". See MiSiPi documentation for descriptions of the weighting methods.
# @param out_type The type of file to write the plots to. Options are "png" or "pdf". Default is PDF.
# @param i The current iteration number
# @param i_total The total number of iterations
# @return plots

.miRNA <- function(chrom_name, reg_start, reg_stop, length, strand,
                   genome_file, bam_file, logfile, wkdir,
                   plot_output, path_to_RNAfold, path_to_RNAplot, write_fastas,
                   weight_reads, out_type, i = NULL, i_total = NULL) {
  
  # i and i_total will be null if called from run_all
  if (!is.null(i)) {
    .inform_iteration(i, i_total, chrom_name, strand)
  }
  
  cat(file = logfile, paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start - 1, " reg_stop: ", reg_stop - 1, "\n"), append = TRUE)
  
  pos <- count <- count.x <- count.y <- end <- r1_end <- r1_start <- dist <- r2_end <- r2_start <- lstop <- lstart <- r1_seq <- loop_seq <- r2_seq <- start <- whole_seq <- width <- NULL

  prefix <- .get_region_string(chrom_name, reg_start, reg_stop)
  
  # do not run locus if length is > 300 - not a miRNA. Also avoids issue where user provides coordinates of miRNA cluster.
  if (reg_stop - reg_start > 300) {
    cat(file = logfile, "length of region is greater than 300. \n", append = TRUE)
    return(.null_mi_res(prefix, strand, wkdir))
  }

  bam_obj <- .open_bam(bam_file, logfile)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)
  chr_name <- names(bam_header[["targets"]])
  chr_length <- unname(bam_header[["targets"]])

  bam_header <- NULL

  # for the read size distribution plot
  #read_dist <- .get_read_dist(bam_obj, chrom_name, reg_start, reg_stop)
  
  stranded_size_dist <- .get_stranded_read_dist(bam_obj, chrom_name, reg_start, reg_stop)
  #.plot_sizes_by_strand(wkdir, stranded_read_dist, chrom_name, reg_start, reg_stop)
  
  chrom <- .get_chr(bam_obj, chrom_name, reg_start, reg_stop, strand)

  ## make the read data tables

  # Processing one strand, so make copy of df to transform
  # 1/9/24 added seq back in for writing miRNA pairs (arms)
  # 7/15/25 Rsamtools::scanBam is returning all reads on the specified strand that overlap the given ScanBamParam region
  # -- This is causing some reads start or end position to be well outside the locus region
  # -- Experimenting with filtering the reads to allow them to be a limited distance beyond the locus region
  
  # Limit reads to 5nt beyond current region of interest
  READ_OVERFLOW_LIMIT <- 5
  
  r2_dt <- .make_si_BamDF(chrom) %>%
    subset(width <= 32 & width >= 18) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1) %>%
    dplyr::filter(
      start >= reg_start - READ_OVERFLOW_LIMIT &
      end <= reg_stop + READ_OVERFLOW_LIMIT) %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n())

  chrom <- NULL

  if (nrow(r2_dt) == 0) {
    return(.null_mi_res(prefix, strand, wkdir))
  } else {
    locus_length <- reg_stop - reg_start + 1
    r2_dt <- .weight_reads(r2_dt, weight_reads, locus_length, sum(r2_dt$count))
  }

  # Transform ends of one set of reads in order capture the other arm of the hairpin
  # Usually no more than 60nt away - Jarva et al lulz
  r1_dt <- r2_dt %>%
    dplyr::mutate(end = end + 59)

  if (nrow(r1_dt) == 0 || nrow(r2_dt) == 0) {
    cat(file = logfile, "After filtering for width and strand, zero reads remain. Please check bam BAM file.\n", append = TRUE)
    return(.null_mi_res(prefix, strand, wkdir))
  }

  pileup_start <- min(r1_dt$start)
  pileup_stop <- max(r2_dt$end)
  pileups <- .get_read_pileups(chrom_name, pileup_start, pileup_stop, strand, bam_file)

  empty_table <- data.frame(pos = c(seq(pileup_start, pileup_stop)), count = c(0))

  dt_table <- merge(empty_table,
    pileups,
    by = "pos",
    all.x = TRUE
  ) %>%
    dplyr::select(-c(count.x)) %>%
    dplyr::rename("count" = count.y)

  dt_table[is.na(dt_table)] <- 0

  #r1_summarized <- r1_dt %>%
  #  dplyr::group_by_all() %>%
  #  dplyr::count()
  
  r1_summarized <- r1_dt %>%
    dplyr::select(-c(first, seq)) %>%
    dplyr::group_by_all() %>%
    dplyr::count()
  
  #r2_summarized <- r2_dt %>%
  #  dplyr::group_by_all() %>%
  #  dplyr::count()
  
  r2_summarized <- r2_dt %>%
    dplyr::select(-c(first, seq)) %>%
    dplyr::group_by_all() %>%
    dplyr::count()

  overlaps <- .find_overlaps(r1_summarized, r2_summarized)

  MAX_DIST_APART <- 60

  overlaps <- overlaps %>%
    dplyr::mutate(r1_end = r1_end - 59) %>%
    dplyr::mutate(r1_width = r1_end - r1_start + 1) %>%
    dplyr::mutate(dist = r2_start - r1_end) %>% # This line and the next effectively replace the cpp function get_nearby()
    dplyr::filter(r1_end < r2_start & dist <= MAX_DIST_APART) %>%
    # dist > 2 to allow for at least a minimal loop
    dplyr::filter(dist > 2) %>% # only want read pairs where start position of read 2 is after the end of read 1
    dplyr::select(-dist)

  # write out pairs here
  if (nrow(overlaps) == 0) {
    cat(file = logfile, "No overlapping reads found.\n", append = TRUE)
    return(.null_mi_res(prefix, strand, wkdir))
  }

  read_pileups <- getPileupsMap(
    dt_table$pos, dt_table$count,
    overlaps$r1_start, overlaps$r1_end,
    overlaps$r2_start, overlaps$r2_end,
    overlaps$r1_dupes, overlaps$r2_dupes
  )

  overlaps <- NULL

  if (nrow(read_pileups) == 0) {
    return(.null_mi_res(prefix, strand, wkdir))
  }

  read_pileups <- read_pileups %>%
    dplyr::mutate(
      width_r1 = (r1_end - r1_start) + 1,
      width_r2 = (r2_end - r2_start) + 1
    )

  read_pileups <- read_pileups %>%
    dplyr::mutate(
      "lstart" = r1_end + 1,
      "lstop" = r2_start - 1
    ) %>% # get rid of rows with loop length less than 3
    dplyr::filter((lstop - lstart + 1) > 2)

  # loop_coord <- loop_coord %>% dplyr::filter((lstop - lstart + 1) > 2)
  if (nrow(read_pileups) == 0) {
    return(.null_mi_res(prefix, strand, wkdir))
  }

  mygranges <- GenomicRanges::GRanges(
    seqnames = c(chrom_name),
    ranges = IRanges::IRanges(start = reg_start - READ_OVERFLOW_LIMIT, end = reg_stop + READ_OVERFLOW_LIMIT)
  )

  bed_seq <- toString(Rsamtools::scanFa(genome_file, mygranges))
  
  # This offset will help give us the relative start and stop positions in the shortened extracted sequence
  RELATIVE_SEQ_POS_OFFSET <- reg_start - READ_OVERFLOW_LIMIT - 1

  loop_seqs <- getFastas(bed_seq, read_pileups$lstart - RELATIVE_SEQ_POS_OFFSET, read_pileups$lstop - RELATIVE_SEQ_POS_OFFSET, nrow(read_pileups)) %>%
    dplyr::rename("lstart" = "start", "lstop" = "stop")

  r1_seqs <- getFastas(bed_seq, read_pileups$r1_start - RELATIVE_SEQ_POS_OFFSET, read_pileups$r1_end - RELATIVE_SEQ_POS_OFFSET, nrow(read_pileups)) %>%
    dplyr::rename("r1_start" = "start", "r1_end" = "stop")
  r2_seqs <- getFastas(bed_seq, read_pileups$r2_start - RELATIVE_SEQ_POS_OFFSET, read_pileups$r2_end - RELATIVE_SEQ_POS_OFFSET, nrow(read_pileups)) %>%
    dplyr::rename("r2_start" = "start", "r2_end" = "stop")


  read_pileups <- read_pileups %>%
    dplyr::mutate(
      w_start = r1_start,
      w_stop = r2_end
    )

  read_pileups$r1_seq <- r1_seqs$Seq
  read_pileups$r2_seq <- r2_seqs$Seq
  read_pileups$loop_seq <- loop_seqs$Seq
  
  loop_seqs <- r1_seqs <- r2_seqs <- NULL
  

  read_pileups$whole_seq <- stringr::str_c(read_pileups$r1_seq, read_pileups$loop_seq, read_pileups$r2_seq)
  read_pileups <- read_pileups %>%
    dplyr::mutate(width = w_stop - w_start + 1) %>%
    dplyr::distinct()


  # select unique combinations of r1_start/stop, r2_start/stop
  grouped <- read_pileups %>%
    dplyr::group_by(r1_start, r1_end, r2_start, r2_end) %>%
    dplyr::distinct() %>%
    dplyr::select(c(r1_start, r1_end, r1_count_avg, r2_start, r2_end, r2_count_avg)) %>%
    dplyr::rename("r1_alt_start" = "r1_start", "r1_alt_end" = "r1_end",
                  "r2_alt_start" = "r2_start", "r2_alt_end" = "r2_end")

  grouped <- grouped %>%
    dplyr::mutate(Chrom = chrom_name, Reg_start = reg_start, Reg_stop = reg_stop) %>%
    dplyr::select(c(Chrom, Reg_start, Reg_stop, r1_alt_start, r1_alt_end, r1_count_avg, r2_alt_start, r2_alt_end, r2_count_avg))

  if (nrow(grouped > 1)) {
    alt_file <- file.path(wkdir, "alt_miRNAs_coord.bed")
    .write.quiet(grouped, alt_file)
    cat(file = logfile, "Writing potential alternative miRNA start and stop coordinates to alt_miRNAs_coord.bed.\n", append = TRUE)
  }
  grouped <- NULL

  most_abundant_idx <- which((read_pileups$r1_count_avg + read_pileups$r2_count_avg) == max(read_pileups$r1_count_avg + read_pileups$r2_count_avg))
  most_abundant <- read_pileups[most_abundant_idx, ]
  
  # Will write this to a file in order to keep track of which strand of each region had the most expression

  # index first value in the even there is more than one most abundant read
  most_abundant_avg_count <- mean(c(most_abundant$r1_count_avg[1], most_abundant$r2_count_avg[1]))

  read_pileups <- NULL
  
  final <- most_abundant[1, ]
  final_seq <- final$whole_seq

  
  roi_seq <- substring(bed_seq, READ_OVERFLOW_LIMIT + 1, nchar(bed_seq) - READ_OVERFLOW_LIMIT)
  # RNAfold wants U's not T's, so convert to U
  expanded_converted <- list(convertU(bed_seq, 1))
  converted <- list(convertU(roi_seq, 1))

  region_string <- paste0(">", chrom_name, "-", reg_start - 1, "_", reg_stop - 1)
  expanded_converted <- data.frame("V1" = unname(unlist(expanded_converted)))
  converted <- data.frame("V1" = unname(unlist(converted)))
  
  # Use bed file coords in column name unless alternate coordinates used
  colnames(expanded_converted) <- region_string
  colnames(converted) <- region_string
  
  # Positions relative to the expanded roi (could use roi bound seq, but the positions would need to be (x - reg_start + 1))
  ma_relative_r1_start <- final$r1_start - RELATIVE_SEQ_POS_OFFSET
  ma_relative_r2_end <- final$r2_end - RELATIVE_SEQ_POS_OFFSET
  most_abundant_seq <- stringr::str_sub(expanded_converted[region_string], ma_relative_r1_start, ma_relative_r2_end)
  final$converted <- most_abundant_seq

  # Use the roi bound converted data frame for converted.fasta for later folding
  write.table(converted, file = file.path(wkdir, "converted.fasta"), sep = "\n", append = FALSE, row.names = FALSE, quote = FALSE)

  mx_idx <- which(c(final$r1_count_avg, final$r2_count_avg) == max(final$r1_count_avg, final$r2_count_avg))
  # mx doesn't appear to be used, and I think the code for this line is wrong
  # If the intention is to create a vector of final$r1_count_avg[mx_idx] and final$r2_count_avg[mx_idx]
  # This will not do that. Instead it should be:
  # mx <- c(final$r1_count_avg[mx_idx], final$r2_count_avg[mx_idx])
  #mx <- c(final$r1_count_avg, final$r2_count_avg)[mx_idx]

  if (length(mx_idx) < 2) {
    # colors <- vector(mode = "character", length = 2)
    colors <- c(NA, NA)
    # set max value color to YELLOW
    colors[mx_idx] <- "GREEN"
    # set not max value to GREEN
    colors[is.na(colors)] <- "RED"
    # colors[which(colors != mx_idx)] <- "GREEN"
  } else {
    colors <- vector(mode = "character", length = 2)
    colors[1] <- "GREEN"
    colors[2] <- "GREEN"
  }

  # need to convert starts and stops of reads to 1-based for coloring structure
  # diff <- final$r1_start - reg_start

  # Get the relative start and stop positions of the reads in the context of final_seq for plotting purposes
  pos_df <- data.frame(
    r1_start = stringr::str_locate(roi_seq, final$r1_seq)[1],
    r1_end = stringr::str_locate(roi_seq, final$r1_seq)[2],
    r2_start = stringr::str_locate(roi_seq, final$r2_seq)[1],
    r2_end = stringr::str_locate(roi_seq, final$r2_seq)[2]
  )

  .rna_plot(path_to_RNAfold, path_to_RNAplot, wkdir, pos_df, colors, chrom_name, reg_start, reg_stop, final$r1_start, final$r2_end, strand)

  ################################################################################################################
  # .fold_short_rna folds a list of sequences whereas fold_long_rna only folds one
  fold_list <- .fold_short_rna(reg_start, reg_stop, converted, path_to_RNAfold, chrom_name, wkdir)
  fold_list$helix <- R4RNA::viennaToHelix(fold_list$vienna)

  # make the plots for all the sequences in the "fold_list"
  mfe <- fold_list$mfe
  perc_paired <- (length(fold_list$helix$i) * 2) / (fold_list$stop - fold_list$start)
  
  # transforms reads from one arm of hairpin to their paired position
  # makes a table of reads which are overlapping
  #dicer_dt, helix_df, chrom_name, reg_start

  dicer_overlaps <- .dicer_overlaps(r2_summarized, fold_list$helix, chrom_name, fold_list$start)
  
  # summarize the counts by the # overlapping nucleotides
  z_res <- make_count_table(r1_dt$start, r1_dt$end, r1_dt$width, r2_dt$start, r2_dt$end, r2_dt$width)
  #z_res <- make_miRNA_count_table(r1_dt$start, r1_dt$end, r1_dt$width, r2_dt$start, r2_dt$end, r2_dt$width)
  
  
  # make_count_table was originally written for piRNAs. Need to subtract 3 from each overlap size.
  # TODO
  # 2/18/25 - NOT SURE AT ALL ABOUT THIS STATEMENT
  # Review the c++ function make_count_table()
  # It is calculating overlaps from 4 through 30
  # I think this is mutation is in correct, and if we want overlaps from 1-30 or 1-27,
  # We'll need to modify the make_count_table to run differently for miRNA than piRNA
  z_res <- z_res %>% dplyr::mutate(overlap = overlap - 3)
  
  # filtered_dcr_overlaps <- dicer_overlaps %>%
  #   dplyr::filter(r1_width >= 18 & r1_width <= 30,
  #                 r2_width >= 18 & r2_width <= 30)
  # 
  # z_res <- make_miRNA_count_table_dcroverlap(
  #   filtered_dcr_overlaps$r2_start,
  #   filtered_dcr_overlaps$r2_end,
  #   filtered_dcr_overlaps$r2_width,
  #   filtered_dcr_overlaps$r2_dupes,
  #   filtered_dcr_overlaps$r1_start,
  #   filtered_dcr_overlaps$r1_end,
  #   filtered_dcr_overlaps$r1_width,
  #   filtered_dcr_overlaps$r1_dupes
  # )
  
  r1_dt <- r2_dt <- NULL
  
  # create empty z_df
  z_df <- data.frame("Overlap" = z_res[, 1], "zscore" = .calc_zscore(z_res$count), "ml_zscore" = .calc_ml_zscore(z_res$count))
  
  z_output <- z_df %>%
    dplyr::select(-ml_zscore)
  
  zdf_output <- as.data.frame(t(z_output[,-1]))
  colnames(zdf_output) <- z_output$Overlap
  
  zdf_file <- file.path(wkdir, "overlap_probability.txt")
  .write.quiet(zdf_output, zdf_file)
  
  # calculate the zscores, if there are results
  if (is.na(dicer_overlaps[1, 1]) | dicer_overlaps[1, 1] == 0) {
    overhangs <- data.frame(shift = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), proper_count = c(0, 0, 0, 0, 0, 0, 0, 0, 0), improper_count = c(0, 0, 0, 0, 0, 0, 0, 0, 0))
    overhangs$zscore <- .calc_zscore(overhangs$proper_count)
    overhangs$ml_zscore <- .calc_ml_zscore(overhangs$proper_count)
  } else {
    # if(write_fastas == TRUE) .write_proper_overhangs(wkdir, prefix, overlaps, "_miRNA")
    overhangs <- calc_overhangs(
      dicer_overlaps$r2_start, dicer_overlaps$r2_end,
      dicer_overlaps$r1_start, dicer_overlaps$r1_width,
      dupes_present = TRUE,
      r1_dupes = dicer_overlaps$r1_dupes,
      r2_dupes = dicer_overlaps$r2_dupes
    )
    
    overhangs$zscore <- .calc_zscore(overhangs$proper_count)
    overhangs$ml_zscore <- .calc_ml_zscore(overhangs$proper_count)
  }

  # transform data frame from table to row
  overhang_output <- data.frame(t(overhangs$ml_zscore))
  colnames(overhang_output) <- overhangs$shift
  overhang_output$original_locus <- prefix
  overhang_output$most_abundant_locus <- .get_region_string(chrom_name, fold_list$start, fold_list$stop)
  overhang_output$strand <- strand
  overhang_output$count_avg <- most_abundant_avg_count
  overhang_output <- overhang_output[, c(10, 11, 12, 13, 1:9)]

  dice_file <- switch(
    strand,
    "+" = "miRNA_plus_dicerz.txt",
    "-" = "miRNA_minus_dicerz.txt"
  )
  dice_file <- file.path(wkdir, dice_file)
  .write.quiet(overhang_output, dice_file)

  if (plot_output == TRUE) {
    plots <- .plot_miRNA(chrom_name, reg_start, reg_stop, strand,
                         bam_file, fold_list, stranded_size_dist,
                         out_type, prefix, wkdir)
  } else {
    plots <- NULL
  }

  results <- list("mfe" = mfe, "perc_paired" = perc_paired, "overhangs" = overhangs, "overlaps" = z_df, "plots" = plots)

  return(results)
}
