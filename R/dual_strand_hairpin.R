# Function to run the hairpin algorithm.
# Processes reads from bam object according to strand.
# Plots the arc plot and read distribution.
# @param chrom_name A string corresponding to the chromosome.
# @param reg_start An integer corresponding to the start of a region of interest.
# @param reg_stop An integer corresponding to the end of a region of interest.
# @param length The length of the chromosome of interest.
# @param genome_file The path to a genome fasta file.
# @param bam_file The path to a BAM file. There must be a corresponding index ending in .bai in the same directory.
# @param logfile The name of the file to which log information will be written.
# @param wkdir The path to the directory where all outputs will be written.
# @param plot_output Determines whether PDF plots will be made. Expected values are TRUE or FALSE
# @param path_to_RNAfold The full path to the RNAfold binary executable.
# @param annotate_region Determines whether the program will plot genomic features of interest found in the GTF annotation file. If TRUE, a GTF file must be provided as the "gtf_file" argument.
# @param weight_reads Determines whether read counts will be weighted and with which method. Valid options are "weight_by_prop", "locus_norm", a user-defined value, or "none". See MiSiPi documentation for descriptions of the weighting methods.
# @param gtf_file A string corresponding to the path of genome annotation in 9-column GTF format.
# @param write_fastas TRUE or FALSE. Default is FALSE
# @param out_type The type of file to write the plots to. Options are "png" or "pdf". Default is PDF.
# @param dicer_overhangs
# @return a list of results

.dual_strand_hairpin <- function(chrom_name, reg_start, reg_stop, length,
                                 genome_file, bam_file, logfile, wkdir, plot_output, path_to_RNAfold, annotate_region,
                                 weight_reads, gtf_file, write_fastas, out_type, dicer_overhangs) {
  end <- dist <- num.y <- num.x <- Zscore <- converted <- NULL


  # function to fold the genomic sequence and extract relevant results
  # calls .fold_long_rna which runs RNAfold from ViennaRNA, splits the strings returned into
  # the dot thing and the mfe
  # R4RNA viennaToHelix is used to create the helix from the dot
  fold_the_rna <- function(geno_seq, chrom_name, reg_start, reg_stop, path_to_RNAfold, wkdir) {
    #dna_vec <- as.character(Biostrings::subseq(geno_seq, start = reg_start, end = reg_stop))
    dna_vec <- as.character(geno_seq)
    converted <- convertU(dna_vec, 1)
    dna_vec <- NULL

    fold_list <- mapply(.fold_long_rna, chrom_name, reg_start, reg_stop, converted, path_to_RNAfold, wkdir)
    
    if(!is.null(unlist(unname(fold_list[1]))) & !is.na(fold_list[4,])){   #for cases where there are no paired bases or other problems from RNAfold
      fold_list <- t(fold_list)
      MFE <- unlist(unname(fold_list[, 3]))
      vienna <- fold_list[, 5]
      extracted_df <- fold_list[4][[1]]
  
      writeLines(as.character(vienna), con = file.path(wkdir, "vienna.txt"))
  
      prefix <- .get_region_string(chrom_name, reg_start, reg_stop)
  
      helix <- R4RNA::viennaToHelix(unlist(fold_list[, 5]))
  
      filepath <- file.path(wkdir, "helix.txt")
      R4RNA::writeHelix(helix, file = filepath)
      return(list(MFE = MFE, vienna = vienna, extracted_df = extracted_df, helix = helix))
    } else {
      res <- .null_hp_res()
      return(res)
    }
  }


  bam_obj <- .open_bam(bam_file, logfile)

  # RNAfold can't fold things longer than 10kb
  if ((reg_stop - reg_start + 1) > 10000) {
    res <- .null_hp_res()
    
    # Add density plot to null result
    data <- .read_densityBySize(chrom_name, reg_start, reg_stop, bam_file, wkdir)
    res$density_plot <- .plot_large_density(data, reg_start, reg_stop)
    
    # Add gtf plot to null result
    res$gtf_plot <- ifelse(annotate_region, .plot_gtf(gtf_file, chrom_name, reg_start, reg_stop, NA))
    
    return(res)
  }

  # Extract chromosome sequence from genome file
  # Bedtools start is 0-indexed, stop is half open 1-based 
  # Need to substract 1 from end to avoid exceeding chromosome length
  mygranges <- GenomicRanges::GRanges(
    seqnames = c(chrom_name),
    ranges = IRanges::IRanges(start = c(reg_start), end = c((reg_stop - 1)))
  )


  geno_seq <- Rsamtools::scanFa(genome_file, mygranges)
  #geno_seq <- as.character(unlist(Biostrings::subseq(geno_seq, start = 1, end = length)))
  geno_seq <- as.character(unlist(geno_seq))
  
  cat(file = logfile, paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start - 1, " reg_stop: ", reg_stop - 1, "\n"), append = TRUE)
  cat(file = logfile, "Filtering forward and reverse reads by length\n", append = TRUE)

  # define which for Rsamtools ScanBamParam
  which <- GenomicRanges::GRanges(seqnames = chrom_name, IRanges::IRanges(reg_start, reg_stop))

  #### Strand Processing ####
  # Process each strand and store results in a shared object for later use
  strands <- c("+", "-")
  sRes <- list()
  sRes$folded <- FALSE
  
  for (strand in strands) {
    null_results <- FALSE
    strand_name <- switch(strand, "+" = "plus", "-" = "minus")
    
    chrom <- .get_chr(bam_obj, chrom_name, reg_start, reg_stop, strand = strand)
    
    r2_dt <- data.table::setDT(.make_si_BamDF(chrom)) %>%
      base::subset(width <= 32 & width >= 18) %>%
      dplyr::rename(start = pos) %>%
      dplyr::mutate(end = start + width - 1) %>%
      dplyr::group_by_all() %>%
      dplyr::summarize(count = dplyr::n())
    
    locus_length <- reg_stop - reg_start + 1
    
    r2_dt <- .weight_reads(r2_dt, weight_reads, locus_length, sum(r2_dt$n))
    
    r2_dt <- na.omit(r2_dt)
    
    # Now that the reads have been weighted,
    # Let's resummarize them for more efficient processing
    r2_dt_summarized <- r2_dt %>%
      dplyr::group_by_all() %>%
      dplyr::count()
    
    # We're operating on a single strand but need one data frame where the reads aren't transformed, one where they are
    # so set the other dt to be the same as the first with transformed ends
    r1_dt_summarized <- r2_dt_summarized %>%
      dplyr::mutate(end = end + 30)
    
    
    # if no results, need to store a specific value for the run_all/machine learning
    # null_hp_res() creates a table of specific "no result" values for zscores and such
    if (nrow(r2_dt) < 3) {
      cat(file = logfile, "After filtering for width and strand, zero reads remain. Please check input BAM file.\n", append = TRUE)
      null_res_idx <- switch(
        strand,
        "-" = 1,
        "+" = 2
      )
      null_res <- .null_hp_res()[[null_res_idx]]
      null_results <- TRUE
    }
    
    # calculate phasing signatures
    if (nrow(r2_dt) == 0) {
      # if read dfs are empty set results to null. Still need to create the empty tables for plots/ML
      cat(file = logfile, "No overlapping reads detected on this strand.\n", append = TRUE)
      
      # creating an empty table with "null" values
      sRes[[strand_name]]$hp_phased_tbl <- data.table::data.table(phased_dist = seq(0, 50), phased_num = rep(0, 51), phased_z = rep(0, 51), phased_ml_z = rep(-33, 51))
      sRes[[strand_name]]$hp_phased_counts <- sum(sRes[[strand_name]]$hp_phased_tbl$phased_num[1:4])
      
      # -33 is an arbitrary value for machine learning purposes
      sRes[[strand_name]]$hp_phased_z <- -33
      sRes[[strand_name]]$hp_phased_mlz <- -33
      sRes[[strand_name]]$perc_paired <- -33
      sRes[[strand_name]]$overhangs <- data.frame(shift = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), proper_count = c(0, 0, 0, 0, 0, 0, 0, 0, 0), improper_count = c(0, 0, 0, 0, 0, 0, 0, 0, 0))
      sRes[[strand_name]]$overhangs$zscore <- .calc_zscore(sRes[[strand_name]]$overhangs$proper_count)
      sRes[[strand_name]]$overhangs$ml_zscore <- .calc_ml_zscore(sRes[[strand_name]]$overhangs$proper_count)
      sRes[[strand_name]]$hp_overhangs_counts <- sum(sRes[[strand_name]]$overhangs$proper_count[5])
      sRes[[strand_name]]$hp_overhangz <- mean(sRes[[strand_name]]$overhangs$zscore[5])
      sRes[[strand_name]]$hp_overhang_mlz <- mean(sRes[[strand_name]]$overhangs$ml_zscore[5])
      sRes[[strand_name]]$all_overlaps <- data.frame()
      
    } else { ## r2_dt has rows
      sRes[[strand_name]]$hp_phased_tbl <- .calc_phasing(r1_dt_summarized, r2_dt_summarized, 50)
      sRes[[strand_name]]$hp_phased_counts <- sum(sRes[[strand_name]]$hp_phased_tbl$phased_num[1:4])
      sRes[[strand_name]]$hp_phased_z <- mean(sRes[[strand_name]]$hp_phased_tbl$phased_z[1:4])
      sRes[[strand_name]]$hp_phased_mlz <- mean(sRes[[strand_name]]$hp_phased_tbl$phased_ml_z[1:4])
      
      # We don't want to fold the dna twice. So once it's been folded
      # set sRes$folded to TRUE
      ## CHECK FOR FOLDED STATUS ##
      if (!sRes$folded) {
        sRes$fold_list <- fold_the_rna(geno_seq, chrom_name, reg_start, reg_stop, path_to_RNAfold, wkdir)
        sRes$MFE <- sRes$fold_list$MFE
        sRes$perc_paired <- (length(sRes$fold_list$helix$i) * 2) / (reg_stop - reg_start)
        sRes$folded <- TRUE
      }
      
      # transform reads and find dicer pairs
      sRes[[strand_name]]$dicer_dt <- r2_dt %>%
        dplyr::group_by(rname, start, end, width, first) %>%
        dplyr::count()
      
      if(!is.null(sRes$fold_list$helix)){ # if RNAfold has problem (e.g. no paired bases, then dicer overlaps function cannot run
        sRes[[strand_name]]$all_overlaps <- .dicer_overlaps(sRes[[strand_name]]$dicer_dt, sRes$fold_list$helix, chrom_name, reg_start)
      
      # .dicer_overlaps() returns zero values if there are no valid overlaps
      # so check to make sure the first values are not zero
      if (!is.na(sRes[[strand_name]]$all_overlaps[1, 1]) && !(nrow(sRes[[strand_name]]$all_overlaps) == 0)) {
        sRes[[strand_name]]$overhangs <- calc_overhangs(
          sRes[[strand_name]]$all_overlaps$r1_start,
          sRes[[strand_name]]$all_overlaps$r1_end,
          sRes[[strand_name]]$all_overlaps$r2_start,
          sRes[[strand_name]]$all_overlaps$r2_width,
          dupes_present = TRUE,
          r1_dupes = sRes[[strand_name]]$all_overlaps$r1_dupes,
          r2_dupes = sRes[[strand_name]]$all_overlaps$r2_dupes
        )
        sRes[[strand_name]]$overhangs$zscore <- .calc_zscore(sRes[[strand_name]]$overhangs$proper_count)
        sRes[[strand_name]]$overhangs$ml_zscore <- .calc_ml_zscore(sRes[[strand_name]]$overhangs$proper_count)
        sRes[[strand_name]]$hp_overhangz <- mean(sRes[[strand_name]]$overhangs$zscore[5])
        sRes[[strand_name]]$hp_overhang_mlz <- mean(sRes[[strand_name]]$overhangs$ml_zscore[5])
      } else {
        sRes[[strand_name]]$overhangs <- data.frame(shift = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), proper_count = c(0, 0, 0, 0, 0, 0, 0, 0, 0), improper_count = c(0, 0, 0, 0, 0, 0, 0, 0, 0))
        sRes[[strand_name]]$overhangs$zscore <- .calc_zscore(sRes[[strand_name]]$overhangs$proper_count)
        sRes[[strand_name]]$overhangs$ml_zscore <- .calc_ml_zscore(sRes[[strand_name]]$overhangs$proper_count)
        # return arbitrary "null" value if there are no valid results for ML
        sRes[[strand_name]]$hp_overhangs_counts <- 0
        sRes[[strand_name]]$hp_overhangz <- 0
        sRes[[strand_name]]$hp_overhang_mlz <- -33
      }
      
      } else {
        sRes[[strand_name]]$all_overlaps <- data.frame()
      }  
        
    }
    
    # results for the ML table
    if (null_results) {
      sRes$MFE <- null_res[[strand_name]]$MFE
      sRes$perc_paired <- null_res[[strand_name]]$perc_paired
      sRes[[strand_name]]$res <- null_res
    } else {
      res <- list(
        MFE = sRes$MFE,
        hp_overhangz = sRes[[strand_name]]$hp_overhangz,
        hp_overhang_mlz = sRes[[strand_name]]$hp_overhang_mlz,
        hp_phasedz = sRes[[strand_name]]$hp_phased_z,
        hp_phased_mlz = sRes[[strand_name]]$hp_phased_mlz,
        phased_tbl.dist = sRes[[strand_name]]$hp_phased_tbl$phased_dist,
        phased_tbl.phased_z = sRes[[strand_name]]$hp_phased_tbl$phased_z,
        phased_tbl.phased_mlz = sRes[[strand_name]]$hp_phased_tbl$phased_ml_z,
        dicer_tbl.shift = sRes[[strand_name]]$overhangs$shift,
        dicer_tbl.zscore = sRes[[strand_name]]$overhangs$zscore,
        dicer_tbl.ml_zscore = sRes[[strand_name]]$overhangs$ml_zscore,
        perc_paired = sRes$perc_paired
      )
      sRes[[strand_name]]$res <- res
    }
  }

  #### Generate Plots ####
  # check and see if dicer results exist for plotting
  if(!nrow(sRes$plus$all_overlaps) == 0){
    plus_overhang_out <- data.frame(t(sRes$plus$res$dicer_tbl.zscore))
    colnames(plus_overhang_out) <- sRes$plus$res$dicer_tbl.shift
    plus_overhang_out$locus <- paste0(chrom_name, "_", reg_start, "_", reg_stop)
    plus_overhang_out <- plus_overhang_out[, c(10, 1:9)]
  
    plus_hp_dicerz_file <- file.path(wkdir, "plus_hp_dicerz.txt")
    .write.quiet(plus_overhang_out, plus_hp_dicerz_file)
  
    minus_overhang_out <- data.frame(t(sRes$minus$res$dicer_tbl.zscore))
    colnames(minus_overhang_out) <- sRes$minus$res$dicer_tbl.shift
    minus_overhang_out$locus <- paste0(chrom_name, "_", reg_start, "_", reg_stop)
    minus_overhang_out <- minus_overhang_out[, c(10, 1:9)]
  
    minus_hp_dicerz_file <- file.path(wkdir, "minus_hp_dicerz.txt")
    .write.quiet(minus_overhang_out, minus_hp_dicerz_file)
  } 
  prefix <- .get_region_string(chrom_name, reg_start, reg_stop)

  plus_phased_out <- t(c(prefix, t(sRes$plus$res$phased_tbl.phased_z)))
  minus_phased_out <- t(c(prefix, t(sRes$minus$res$phased_tbl.phased_z)))

  plus_hp_phasedz_file <- file.path(wkdir, "plus_hp_phasedz.txt")
  minus_hp_phasedz_file <- file.path(wkdir, "minus_hp_phasedz.txt")

  .write.quiet(plus_phased_out, plus_hp_phasedz_file)
  .write.quiet(minus_phased_out, minus_hp_phasedz_file)

  # Base results for machine learning that get returned regardless of plot_output status
  results <- list(
    minus_res = sRes$minus$res,
    plus_res = sRes$plus$res
  )

  # Add additional data and plots to results
  if (plot_output) {
    
    if(!nrow(sRes$plus$all_overlaps) == 0) {
      plus_overhangs <- data.frame(shift = sRes$plus$res$dicer_tbl.shift, zscore = sRes$plus$res$dicer_tbl.zscore)
      minus_overhangs <- data.frame(shift = sRes$minus$res$dicer_tbl.shift, zscore = sRes$minus$res$dicer_tbl.zscore)
  
      ## return these as plot objects
      #plus_overhang_plot <- .plot_overhangz(plus_overhangs, "+")
      #minus_overhang_plot <- .plot_overhangz(minus_overhangs, "-")
      overhang_probability_plot <- .plot_siRNA_overhangs_combined(plus_overhangs, minus_overhangs, dicer_overhangs)
    } 
    
    data <- .readDensityBySize(chrom_name, reg_start, reg_stop, bam_file, wkdir)

    density_plot <- .plot_density(data, reg_start, reg_stop)

    # Check to see if helix.txt has actual data
    helix_file <- file.path(wkdir, "helix.txt")
    
        
    if (file.exists(helix_file)) {
      # First line in helix.txt is just a comment, like # 22, so skip this line
      invisible(nrows_helix_file <- nrow(readr::read_tsv(helix_file, skip = 1, show_col_types = FALSE)))

      if(nrows_helix_file > 0){
        if (!Sys.info()["sysname"] == "Windows") {
          arc_plot <- .plot_helix(helix_file)
          grDevices::dev.control("enable")
          R4RNA::plotHelix(helix = R4RNA::readHelix(helix_file), line = TRUE, arrow = FALSE, lwd = 2.25, scale = FALSE)
          arc_plot <- grDevices::recordPlot() # don't touch this...the boss gets mad
        } else {
          # Read helix data
          helix_data <- R4RNA::readHelix(helix_file)
          
          plot_arc2 <- function() {
            R4RNA::plotHelix(
              helix = helix_data,
              line = TRUE,
              arrow = FALSE,
              lwd = 2.25,
              scale = FALSE
            )
            title(main = "RNAfold Arc", line = -3, font.main = 1)
          }
          
          arc_plot <- gridGraphics::echoGrob(plot_arc2)
          
        }
      }
    } else {
      arc_plot <- NA
    }
    

    ## why? No one knows

    phasedz <- .plot_siRNA_hp_phasing_probability_combined(sRes$plus$hp_phased_tbl, sRes$minus$hp_phased_tbl)
    #plus_phasedz <- .plot_hp_phasedz(sRes$plus$hp_phased_tbl, "+")
    #minus_phasedz <- .plot_hp_phasedz(sRes$minus$hp_phased_tbl, "-")

    if(!nrow(sRes$plus$all_overlaps) == 0){
      results$overhang_probability_plot <- overhang_probability_plot
      #results$plus_overhang_plot = plus_overhang_plot
      #results$minus_overhang_plot = minus_overhang_plot
    } else {
      results$overhang_probability_plot <- NA
      #results$plus_overhang_plot <- NA
      #results$minus_overhang_plot <- NA
    }
    results$density_plot <- density_plot
    results$arc_plot <- arc_plot
    results$phasedz <- phasedz

    # Plot genome annotations (optional)
    results$gtf_plot <- ifelse(annotate_region, .plot_gtf(gtf_file, chrom_name, reg_start, reg_stop), NA)

  return(results)
  }
}
