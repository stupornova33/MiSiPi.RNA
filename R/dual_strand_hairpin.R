dual_strand_hairpin <- function(
    chrom_name, reg_start, reg_stop, genome_file, prefix, locus_length, 
    dicer_df_plus, dicer_df_minus, f_df, r_df, path_to_RNAfold, logfile, wkdir, density_timeout) {
  
  # Return null results if locus too long to fold
  if (locus_length > 10000) {
    output_msg <- paste0("Large locus size: ", locus_length, " for locus: ", prefix, ". Unable to fold.\n")
    cat(file = logfile, output_msg, append = TRUE)
    
    fold_list <- NA
    MFE <- -33
    perc_paired <- 0
    plus_dsh <- .null_hp_res()
    minus_dsh <- .null_hp_res()
    combined_dsh <- list(
      plus_dsh = plus_dsh,
      minus_dsh = minus_dsh,
      MFE = MFE, 
      perc_paired = perc_paired)
    return(combined_dsh)
  }
  
  # Extract region sequence from genome
  mygranges <- GenomicRanges::GRanges(
    seqnames = c(chrom_name),
    ranges = IRanges::IRanges(start = reg_start, end = reg_stop)
  )
  
  geno_seq <- Rsamtools::scanFa(genome_file, mygranges)
  geno_seq <- as.character(unlist(geno_seq))
  
  fold_list <- fold_rna(geno_seq, prefix, path_to_RNAfold, wkdir)
  geno_seq <- NULL
  
  if (fold_list$results_present == FALSE) {
    # Return null results
    plus_dsh <- .null_hp_res()
    minus_dsh <- .null_hp_res()
    combined_dsh <- list(
      plus_dsh = plus_dsh,
      minus_dsh = minus_dsh,
      MFE = MFE,
      perc_paired = 0)
    
    return(combined_dsh)
  }
  
  # Results present from folding
  perc_paired <- (length(fold_list$helix$i) * 2) / locus_length
  MFE <- fold_list$MFE
  
  if (nrow(f_df) == 0) {
    output_msg <- "No reads available on forward strand. Setting null hairpin results.\n"
    cat(file = logfile, output_msg, append = TRUE)
    plus_dsh <- .null_hp_res()
  } else {
    plus_dsh <- process_hairpin_strand(chrom_name, reg_start, dicer_df_plus, f_df, fold_list, logfile, wkdir, density_timeout)
  }
  
  if (nrow(r_df) == 0) {
    output_msg <- "No reads available on reverse strand. Setting null hairpin results.\n"
    cat(file = logfile, output_msg, append = TRUE)
    minus_dsh <- .null_hp_res()
  } else {
    minus_dsh <- process_hairpin_strand(chrom_name, reg_start, dicer_df_minus, r_df, fold_list, logfile, wkdir, density_timeout)
  }
  
  # write results to files
  output_dsh(plus_dsh, minus_dsh, prefix, wkdir)
  
  combined_dsh <- list(
    plus_dsh = plus_dsh,
    minus_dsh = minus_dsh,
    MFE = MFE,
    perc_paired = perc_paired)
  
  return(combined_dsh)
}

# Function to run the hairpin algorithm.
# Processes reads from bam object according to strand.
# @param chrom_name A string corresponding to the chromosome.
# @param reg_start An integer corresponding to the start of a region of interest.
# @param dicer_df A data frame of grouped and summarized reads
# @param df_summarized A data frame of weighted, grouped and summarized reads
# @param fold_list The return object from Vienna RNAfold
# @param logfile The name of the file to which log information will be written.
# @param wkdir The path to the directory where all outputs will be written.
# @param density_timeout The amount of time to try to calculate phasing before returning null results. Default 3600 
# @return a list of results

# New dsh function
process_hairpin_strand <- function(chrom_name, reg_start, dicer_df, df_summarized, fold_list, logfile, wkdir, density_timeout) { 
  r2_dt <- df_summarized
  # We're operating on a single strand but need one data frame where the reads aren't transformed, one where they are
  # so set the other dt to be the same as the first with transformed ends
  r1_dt <- r2_dt %>%
    dplyr::mutate(end = end + 30)
  
  hp_phased_tbl <- NULL 
  
  timeout <- density_timeout 
  
  tryCatch( 
    #### Phasing Signatures #### 
    # should n being passed in be 30, not 50? 
    hp_phased_tbl <- R.utils::withTimeout( 
      .calc_phasing(r1_dt, r2_dt, 50), 
      # hp_phased_counts <- sum(hp_phased_tbl$phased_num[1:4]) # Appears to not be used anywhere 
      timeout = timeout 
    ),  
    # Log the errors and move on  
    error = function(e) {  
      error_msg <- paste("ERROR: Timeout of", timeout, "seconds was exceeded by .calc_phasing from dual_strnd_hairpin. Timeout can be increased in set_vars.\nSetting null results for phasing and continuing.")  
      cat(file = logfile, error_msg)  
    }
  ) 
  
  if (!is.null(hp_phased_tbl)) { 
    hp_phased_z <- mean(hp_phased_tbl$phased_z[1:4]) 
    hp_phased_mlz <- mean(hp_phased_tbl$phased_ml_z[1:4]) 
  } else { 
    tmp <- .null_hp_res() 
    hp_phased_tbl <- tmp$hp_phased_tbl 
    hp_phased_z <- tmp$hp_phasedz 
    hp_phased_mlz <- tmp$hp_phased_mlz 
  } 
  
  r1_dt <- r2_dt <- NULL
  
  #### Dicer Overlaps ####
  if (is.null(fold_list$helix)) {
    # NULL all_overlaps and overhangs when helix not present
    all_overlaps <- data.frame()
    overhangs <- data.frame(
      shift = -4:4,
      proper_count = 0,
      improper_count = 0,
      zscore = 0,
      ml_zscore = -33,
      hp_overhangz = 0,
      hp_overhang_mlz <- -33
    )
    
  } else {
    # Helix present so calculate dicer overlaps
    all_overlaps <- .dicer_overlaps(dicer_df, fold_list$helix, reg_start)
    
    # Check for 0 values indicating no valid overlaps
    if (is.na(all_overlaps[1, 1]) | nrow(all_overlaps) == 0) {
      # NULL overhangs when no valid overlaps
      overhangs <- data.frame(
        shift = -4:4,
        proper_count = 0,
        improper_count = 0,
        zscore = 0,
        ml_zscore = -33,
        hp_overhangz = 0,
        hp_overhang_mlz <- -33
      )
      
    # Valid overlaps are present so calculate overhangs
    } else {
      overhangs <- calc_overhangs(
        all_overlaps$r1_start,
        all_overlaps$r1_end,
        all_overlaps$r2_start,
        all_overlaps$r2_width,
        dupes_present = TRUE,
        r1_dupes = all_overlaps$r1_dupes,
        r2_dupes = all_overlaps$r2_dupes
      )
      overhangs$zscore <- .calc_zscore(overhangs$proper_count)
      overhangs$ml_zscore <- .calc_ml_zscore(overhangs$proper_count)
      overhangs$hp_overhangz <- mean(overhangs$zscore[5])
      overhangs$hp_overhang_mlz <- mean(overhangs$ml_zscore[5])
    }
  }
  # TODO Add in calc zscore for overhangs: see line 190 in dual_strand_haipin
  results <- list(
    all_overlaps = all_overlaps,
    hp_overhangz = overhangs$hp_overhangz,
    hp_overhang_mlz = overhangs$hp_overhang_mlz,
    hp_phasedz = hp_phased_z, # Avg of hp_phased_tbl first 4 z scores
    hp_phased_mlz = hp_phased_mlz, # Avg of hp_phased_tbl first 4 ml z scores
    hp_phased_tbl = hp_phased_tbl,
    phased_tbl.dist = hp_phased_tbl$phased_dist,
    phased_tbl.phased_z = hp_phased_tbl$phased_z,
    phased_tbl.phased_mlz = hp_phased_tbl$phased_ml_z,
    dicer_tbl.shift = overhangs$shift,
    dicer_tbl.zscore = overhangs$zscore,
    dicer_tbl.ml_zscore = overhangs$ml_zscore
  )

  return(results)
}

output_dsh <- function(dsh_plus, dsh_minus, prefix, wkdir) {
  if (nrow(dsh_plus$all_overlaps) != 0) {
    #### Plus strand hairpin dicer zscores output file ####
    plus_overhang_out <- data.frame(t(dsh_plus$dicer_tbl.zscore))
    colnames(plus_overhang_out) <- dsh_plus$dicer_tbl.shift
    plus_overhang_out$locus <- prefix
    plus_overhang_out <- plus_overhang_out[, c(10, 1:9)]
    
    plus_hp_dicerz_file <- file.path(wkdir, "plus_hp_dicerz.txt")
    .write.quiet(plus_overhang_out, plus_hp_dicerz_file)
  }
  
  if (nrow(dsh_minus$all_overlaps) != 0) {
    #### Minus strand hairpin dicer zscores output file ####
    minus_overhang_out <- data.frame(t(dsh_minus$dicer_tbl.zscore))
    colnames(minus_overhang_out) <- dsh_minus$dicer_tbl.shift
    minus_overhang_out$locus <- prefix
    minus_overhang_out <- minus_overhang_out[, c(10, 1:9)]
    
    minus_hp_dicerz_file <- file.path(wkdir, "minus_hp_dicerz.txt")
    .write.quiet(minus_overhang_out, minus_hp_dicerz_file)
  }
  
  #### Plus strand hairpin phased zscores output file ####
  plus_phased_out <- t(c(prefix, t(dsh_plus$phased_tbl.phased_z)))
  plus_hp_phasedz_file <- file.path(wkdir, "plus_hp_phasedz.txt")
  .write.quiet(plus_phased_out, plus_hp_phasedz_file)
  
  #### Minus strand hairpin phased zscores output file ####
  minus_phased_out <- t(c(prefix, t(dsh_minus$phased_tbl.phased_z)))
  minus_hp_phasedz_file <- file.path(wkdir, "minus_hp_phasedz.txt")
  .write.quiet(minus_phased_out, minus_hp_phasedz_file)
}

plot_dsh <- function(dsh_plus, dsh_minus, dicer_overhangs, wkdir) {
  #### Overhang Probability Plot ####
  if (nrow(dsh_plus$all_overlaps) == 0) {
    overhang_probability_plot <- null_plot("siRNA Dicer Overhang Probability", "No overlaps were present")
    
  } else {
    plus_overhangs <- data.frame(shift = dsh_plus$dicer_tbl.shift, zscore = dsh_plus$dicer_tbl.zscore)
    minus_overhangs <- data.frame(shift = dsh_minus$dicer_tbl.shift, zscore = dsh_minus$dicer_tbl.zscore)
    
    overhang_probability_plot <- .plot_siRNA_overhangs_combined(plus_overhangs, minus_overhangs, dicer_overhangs)
  }
  
  #### Arc Plot ####
  helix_file <- file.path(wkdir, "helix.txt")
  
  if (!file.exists(helix_file)) {
    arc_plot <- null_plot("RNAfold Arc Plot", "Not generated due to RNA not being folded.")
  } else {
    invisible(nrows_helix_file <- nrow(readr::read_tsv(helix_file, skip = 1, show_col_types = FALSE)))
    
    if (nrows_helix_file == 0) {
      arc_plot <- null_plot("RNAfold Arc Plot", "Not generated due to empty helix file.")
    } else {
      arc_plot <- plot_helix(helix_file)
    }
  }
  
  #### Hairpin Phasing Probability Plot ####
  if (dsh_plus$hp_phasedz == 0 & dsh_minus$hp_phasedz == 0) {
    phasedz <- null_plot("siRNA Phasing Probability", "No results on which to test phasing")
  } else {
    phasedz <- .plot_siRNA_hp_phasing_probability_combined(dsh_plus$hp_phased_tbl, dsh_minus$hp_phased_tbl)
  }
  
  #### Gather and return plots ####
  plots <- list(
    overhang_probability_plot = overhang_probability_plot,
    arc_plot = arc_plot,
    phasedz = phasedz)
  
  return(plots)
}

plot_helix <- function(helix_file) {
  # Read in Helix data
  helix_data <- R4RNA::readHelix(helix_file)
  
  plot_arc <- function() {
    R4RNA::plotHelix(
      helix = helix_data,
      line = TRUE,
      arrow = FALSE,
      lwd = 2.25,
      scale = FALSE
    )
    title(main = "RNAfold Arc", line = -3, font.main = 1)
  }
  
  # Plot generation specific to user operating system
  platform <- Sys.info()["sysname"]
  
  if (platform == "Windows") {
    arc_plot <- gridGraphics::echoGrob(plot_arc)
    
  } else { # Linux / Unix based
    grDevices::dev.control("enable")
    plot_arc()
    arc_plot <- grDevices::recordPlot()
  }
  
  return(arc_plot)
}

# Fold the genomic sequence and extract relevant results
# calls .fold_long_rna which runs RNAfold from ViennaRNA, splits the strings returned into
# the dot thing and the mfe
# R4RNA viennaToHelix is used to create the helix from the dot
fold_rna <- function(geno_seq, prefix, path_to_RNAfold, wkdir) {
  converted <- convertU(geno_seq, 1)
  
  fold_list <- .fold_long_rna(prefix, converted, path_to_RNAfold, wkdir)
  
  # new fold_long_rna return structure
  # [1] - mfe
  # [2] - vien_struct
  
  MFE <- fold_list$mfe
  vienna <- fold_list$vien_struct
  helix <- NA
  
  # Default value for checking return status later
  results_present <- FALSE
  
  # Create helix if results present
  # NA is default value for when results are not present
  if (!is.na(vienna)) {
    results_present <- TRUE
    helix <- R4RNA::viennaToHelix(vienna)
    
    writeLines(as.character(vienna), con = file.path(wkdir, "vienna.txt"))
    
    helix_filepath <- file.path(wkdir, "helix.txt")
    R4RNA::writeHelix(helix, file = helix_filepath)
  }
  
  result <- list(
    results_present = results_present,
    MFE = MFE,
    vienna = vienna,
    helix = helix)
  
  return(result)
}
