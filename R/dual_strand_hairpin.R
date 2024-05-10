#' Function to run the hairpin algorithm.
#' Processes reads from bam object according to strand.
#' Plots the arc plot and read distribution.
#' @param chrom_name A string corresponding to the chromosome.
#' @param reg_start An integer corresponding to the start of a region of interest.
#' @param reg_stop An integer corresponding to the end of a region of interest.
#' @param length The length of the chromosome of interest.
#' @param min_read_count A whole number. Default is 1.
#' @param genome_file The path to a genome fasta file.
#' @param bam_file The path to a BAM file. There must be a corresponding index ending in .bai in the same directory.
#' @param logfile The name of the file to which log information will be written.
#' @param wkdir The path to the directory where all outputs will be written.
#' @param plot_output Determines whether PDF plots will be made. Expected values are "T" or "F".
#' @param path_to_RNAfold The full path to the RNAfold binary executable.
#' @param annotate_region Determines whether the program will plot genomic features of interest found in the GTF annotation file. If "T", a GTF file must be provided as the "gtf_file" argument.
#' @param weight_reads Determines whether read counts will be weighted and with which method. Valid options are "weight_by_prop", "locus_norm", a user-defined value, or "none". See MiSiPi documentation for descriptions of the weighting methods.
#' @param gtf_file A string corresponding to the path of genome annotation in 9-column GTF format.
#' @param out_type The type of file to write the plots to. Options are "png" or "pdf". Default is PDF.
#' @return a list of results

#' @export

dual_strand_hairpin <- function(chrom_name, reg_start, reg_stop, length,
                             min_read_count, genome_file, bam_file, logfile, wkdir, plot_output, path_to_RNAfold, annotate_region,
                             weight_reads, gtf_file, out_type){

  end <- dist <- num.y <- num.x <- Zscore <- converted <- NULL


  # function to fold the genomic sequence and extract relevant results
  # calls fold_long_rna which runs RNAfold from ViennaRNA, splits the strings returned into
  # the dot thing and the mfe
  # R4RNA viennaToHelix is used to create the helix from the dot
  fold_the_rna <- function(geno_seq, chrom_name, reg_start, reg_stop, path_to_RNAfold){

     dna_vec <- as.character(Biostrings::subseq(geno_seq, start = reg_start, end = reg_stop))

     converted <- convertU(dna_vec, 1)
     writeLines(converted, con = "converted.txt")
     dna_vec <- NULL

     # use
     fold_list <- mapply(fold_long_rna, chrom_name, reg_start, reg_stop, converted, path_to_RNAfold)
     fold_list <- t(fold_list)
     MFE <- unlist(unname(fold_list[,3]))
     vienna <- fold_list[,5]
     extracted_df <- fold_list[4][[1]]

     writeLines(as.character(vienna), con = "vienna.txt")

     prefix <- paste0(wkdir, chrom_name, "-", reg_start, "_", reg_stop, "_", strand)

     helix <- R4RNA::viennaToHelix(unlist(fold_list[,5]))
     R4RNA::writeHelix(helix, file = "helix.txt")
     return(list(MFE = MFE, vienna = vienna, extracted_df = extracted_df, helix = helix))
  }


  bam_obj <- OpenBamFile(bam_file, logfile)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)

  #RNAfold can't fold things longer than 10kb

  if(reg_stop - reg_start > 10000){
    print("Region greater than 10kb. Creating null_hp_res.")
    res <- null_hp_res()
    return(res)
  }
  bam_header <- NULL

  # Extract chromosome sequence from genome file
  mygranges <- GenomicRanges::GRanges(
    seqnames = c(chrom_name),
    ranges = IRanges::IRanges(start=c(1), end=c(length)))

  #mygranges <- GenomicRanges::GRanges(seqnames = c(chrom_name), ranges = IRanges::IRanges(start = c(reg_start),
  #                                                                                        end = c(reg_stop)))
  geno_seq <- Rsamtools::scanFa(genome_file, mygranges)
  geno_seq <- as.character(unlist(Biostrings::subseq(geno_seq, start = 1, end = length)))

  cat(file = paste0(wkdir, logfile), paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start, " reg_stop: ", reg_stop, "\n"), append = TRUE)
  cat(file = paste0(wkdir, logfile), "Filtering forward and reverse reads by length\n", append = TRUE)

  #define which for Rsamtools ScanBamParam
  which <- GenomicRanges::GRanges(seqnames=chrom_name, IRanges::IRanges(reg_start, reg_stop))

  ############################################################ compute plus strand ########################################################
  print("Computing plus strand.")
  strand <- "+"
  bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = FALSE), what=c('rname', 'pos', 'qwidth'), which=which)
  chrom <- getChrPlus(bam_obj, chrom_name, reg_start, reg_stop)


  print("Filtering reads")
  #filter the reads and calculate the end position
  filter_r2_dt <- data.table::setDT(make_si_BamDF(chrom)) %>%
    base::subset(width <= 32 & width >= 18) %>%
    dplyr::mutate(start = pos, end = pos + width - 1) %>%
    dplyr::select(-c(pos)) %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n())

  # We're operating on a single strand but need one data frame where the reads aren't transformed, one where they are
  # so set the other dt to be the same as the first
  dt <- filter_r2_dt

  print("Weighting reads.")
  # weight reads if argument supplied
  if(weight_reads == "weight_by_prop"){
    r2_dt <- weight_by_prop(dt, chrom_name)
  } else if(weight_reads == "Locus_norm" | weight_reads == "locus_norm"){
    r2_dt <- locus_norm(dt, sum(dt$count))
  } else if(is.integer(weight_reads)){
    r2_dt <- weight_by_uservalue(dt, weight_reads, (reg_stop - reg_start)) %>% dplyr::mutate(width = end - start + 1)
  } else {
    r2_dt <- no_weight(dt, chrom_name)
  }

  #transform end of reads in one df
  print("Transforming ends of reads.")
  r1_dt <- r2_dt %>% dplyr::mutate(end = end + 30)

  # if no results, need to store a specific value for the run_all/machine learning
  # null_hp_res() creates a table of specific "no result" values for zscores and such
  print("r1_dt contains less than 3 rows. Setting plus_null_res.")
  if(nrow(r1_dt) < 3 || nrow(r2_dt) < 3){
    cat(file = paste0(wkdir, logfile), "After filtering for width and strand, zero reads remain. Please check input BAM file.\n", append = TRUE)
    plus_null_res <- null_hp_res()[[2]]
  }

  # calculate phasing signatures
  print("Calculating plus phasing signature.")
  if(nrow(r1_dt) > 0 && nrow(r2_dt) > 0){
    plus_hp_phased_tbl <- calc_phasing(r1_dt, r2_dt, 30)
    plus_hp_phased_counts <- sum(plus_hp_phased_tbl$phased_num[1:4])
    plus_hp_phased_z <- mean(plus_hp_phased_tbl$phased_z[1:4])
  } else {
    # if read dfs are empty set results to null. Still need to create the empty tables for plots/ML
    cat(file = paste0(wkdir, logfile), "No overlapping reads detected on this strand.\n", append = TRUE)
    # creating an empty table with "null" values
    plus_hp_phased_tbl <- data.table::data.table(phased_dist = seq(0,50), phased_num = rep(0,51), zscore = rep(0,51))
    plus_hp_phased_counts <- sum(plus_hp_phased_tbl$phased_num[1:4])
    # -33 is an arbitrary value
    plus_hp_phased_z <- -33
  }

  print("Checking to see if result is NA.")
  if(plus_hp_phased_z == "NaN"){
    plus_hp_phased_z <- -33
  }


  if(nrow(r2_dt) > 0){
    print("r2_dt contains data. Proceeding with fold and overlap calc.")
    # don't want to fold the dna twice. So once it's been folded
    # set fold_bool to TRUE
    fold_bool <- 'TRUE'
    print("Folding the RNA.")
    fold_list <- fold_the_rna(geno_seq, chrom_name, reg_start, reg_stop, path_to_RNAfold)
    MFE <- fold_list$MFE
    perc_paired <- (length(fold_list$helix$i)*2)/(reg_stop - reg_start)

    #transform reads and find dicer pairs
    print("Calculating dicer overlaps.")
    #system.time(all_overlaps <- dicer_overlaps(r2_dt, fold_list$helix, chrom_name, reg_start))
    r2_dt <- r2_dt %>% dplyr::group_by(rname, start, end, width, first) %>% dplyr::summarize(count = dplyr::n())

    all_overlaps <- dicer_overlaps(r2_dt, fold_list$helix, chrom_name, reg_start)
    # dicer_overlaps() returns zero values if there are no valid overlaps
    # so check to make sure the first values are not zero
    if(!is.na(all_overlaps[1,1]) && !(all_overlaps[1,1] == 0)){
      print("all_overlaps contains data. Calculating overhangs.")
      plus_overhangs <- data.frame(calc_overhangs(all_overlaps$r1_start, all_overlaps$r1_end,
                                                  all_overlaps$r2_start, all_overlaps$r2_width))
      plus_overhangs$zscore <- calc_zscore(plus_overhangs$proper_count)
      plus_hp_overhangz <- mean(plus_overhangs$zscore[5])

    } else {
      print("all_overlaps does not contain data. Setting null results.")
      plus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
      plus_overhangs$zscore <- calc_zscore(plus_overhangs$proper_count)
      # return arbitrary "null" value if there are no valid results for ML
      plus_hp_overhangs_counts <- 0
      plus_hp_overhangz <- -33
    }
  } else {
      print("r2_dt contains less than 2 rows. Setting null results.")
      plus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
      plus_overhangs$zscore <- calc_zscore(plus_overhangs$proper_count)
      plus_hp_overhangs_counts <- sum(plus_overhangs$proper_count[5])
      plus_hp_overhangz <- mean(plus_overhangs$zscore[5])

      # if there were no results, and the dna didn't get folded, set fold_bool to FALSE so it gets folded in the minus strand part
      fold_bool <- 'FALSE'
      perc_paired <- -33
  }


  # results for the ML table
  print("Creating plus_res")
  if(exists("plus_null_res")){
    print("plus_null_res exists.")
    plus_res <- plus_null_res
  } else {
    print("plus_null_res does not exist. Creating final plus_res object.")
      plus_overhangs$zscore <- calc_zscore(plus_overhangs$proper_count)
      plus_overhangz <- mean(plus_overhangs$zscore[1:4])
      plus_res <- list(plusMFE = MFE, plus_hp_overhangz = plus_hp_overhangz, plus_hp_phasedz = plus_hp_phased_z, phased_tbl.dist = plus_hp_phased_tbl$phased_dist,
                phased_tbl.zscore = plus_hp_phased_tbl$phased_z, dicer_tbl.shift = plus_overhangs$shift, dicer_tbl.zscore = plus_overhangs$zscore, perc_paired= perc_paired)
  }
  ############################################################# compute minus strand ############################################################
  # do the same thing for the minus strand
  print("Beginning minus strand.")
  strand <- "-"
  bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = TRUE), what=c('rname', 'pos', 'qwidth'), which=which)
  chrom <- getChrMinus(bam_obj, chrom_name, reg_start, reg_stop)

  print("Filtering data.")
  filter_r2_dt <- data.table::setDT(make_si_BamDF(chrom)) %>%
    base::subset(width <= 32 & width >= 18) %>%
    dplyr::mutate(start = pos, end = pos + width - 1) %>%
    dplyr::select(-c(pos)) %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n())

  dt <- filter_r2_dt
  print("Weighting reads.")
  # weight reads if argument supplied
  if(weight_reads == "weight_by_prop"){
    r2_dt <- weight_by_prop(dt, chrom_name)
  } else if(weight_reads == "Locus_norm" | weight_reads == "locus_norm"){
    r2_dt <- locus_norm(dt, sum(dt$count))
  } else if(is.integer(weight_reads)){
    r2_dt <- weight_by_uservalue(dt, weight_reads, (reg_stop - reg_start)) %>% dplyr::mutate(width = end - start + 1)
  } else {
    r2_dt <- no_weight(dt, chrom_name)
  }
  #transform ends of reads for phasing/finding overlaps
  print("Transforming ends of reads.")
  r1_dt <- r2_dt %>% dplyr::mutate(end = end + 30)

  if(nrow(r1_dt) < 3 || nrow(r2_dt) < 3){
    print("nrow r1_dt is less than 3. Setting null_minus_res.")
    cat(file = paste0(wkdir, logfile), "After filtering for width and strand, zero reads remain. Please check input BAM file.\n", append = TRUE)
    minus_null_res <- null_hp_res()[[1]]

  }

  #calculate phasing
  if(nrow(r1_dt) > 0 && nrow(r2_dt) > 0){
     print("r1_dt contains data. Calculating phasing.")
     minus_hp_phased_tbl <- calc_phasing(r1_dt, r2_dt, 30)
     print("summing minus_hp phased_num.")
     minus_hp_phased_counts <- sum(minus_hp_phased_tbl$phased_num[1:4])
     print("getting mean of minus_hp phased_num.")
     minus_hp_phasedz <- mean(minus_hp_phased_tbl$phased_z[1:4])
     print("finished getting mean of phased_num.")
  } else {
    print("r1_dt does not contain data. Setting null phasing results.")
    cat(file = paste0(wkdir, logfile), "No overlapping reads detected on this strand.\n", append = TRUE)
    minus_hp_phased_tbl <- data.table::data.table(phased_dist = seq(0,50), phased_num = rep(0,51), phased_z = rep(0,51))
    minus_hp_phased_counts <- sum(minus_hp_phased_tbl$phased_num[1:4])
    minus_hp_phasedz <- -33
  }

  print("minus_hp_phasedz is NA. Setting to -33.")
  if(minus_hp_phasedz == "NaN"){
    minus_hp_phasedz <- -33
  }

  #i 9 j 1
  print("Beginning dicer stuff.")
  if(nrow(r2_dt) > 0){
    print("nrow r2_dt > 0.")
    print(paste0("fold_bool: ", fold_bool))
    if(fold_bool == 'TRUE'){
       print("Calculating dicer_overlaps.")

       # take unique reads for dicer overhang calculation, then replicate according to count later
       r2_dt <- r2_dt %>% dplyr::group_by(rname, start, end,first, width) %>% dplyr::summarize(count = dplyr::n())
       #system.time(all_overlaps <- dicer_overlaps(r2_dt, fold_list$helix, chrom_name, reg_start))
       all_overlaps <- dicer_overlaps(r2_dt, fold_list$helix, chrom_name, reg_start)

       if(!is.na(all_overlaps[1,1]) && !(all_overlaps[1,1] == 0)){  #if there are overlaps calc overhangs
         print("all_overlaps contains results. Computing overhangs.")
         minus_overhangs <- data.frame(calc_overhangs(all_overlaps$r1_start, all_overlaps$r1_end,
                                                     all_overlaps$r2_start, all_overlaps$r2_width))
         minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
         minus_hp_overhangz <- mean(plus_overhangs$zscore[5])
       } else {
           print("all_overlaps does not contain results. Setting results to null.")
           minus_hp_overhangs_counts <- 0
           minus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
           minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
           minus_hp_overhangs_counts <- sum(minus_overhangs$proper_count[5])
           minus_hp_overhangz <- mean(minus_overhangs$zscore[5])
       }
    } else { #else if fold bool is false and r2_dt > 0
        print("Fold_bool == FALSE.")
        fold_list <- fold_the_rna(geno_seq, chrom_name, reg_start, reg_stop, path_to_RNAfold)
        MFE <- fold_list$MFE
        perc_paired <- (length(fold_list$helix$i)*2)/(reg_stop - reg_start)
        all_overlaps <- dicer_overlaps(r2_dt, fold_list$helix, chrom_name, reg_start)

      if(!is.na(all_overlaps[1,1]) && !(all_overlaps[1,1] == 0)){  #if there are overlaps calc overhangs
         minus_overhangs <- data.frame(calc_overhangs(all_overlaps$r1_start, all_overlaps$r1_end,
                                                      all_overlaps$r2_start, all_overlaps$r2_width))
         minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
         minus_hp_overhangz <- mean(minus_overhangs$zscore[5])
      } else {
         minus_hp_overhangs_counts <- 0
         minus_hp_overhangz <- -33

         minus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
         minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
      }
   }

  } else { #else if fold bool is false and no results in r2_dt
      print("fold_bool == FALSE & nrow r2_dt == 0.")
      minus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
      minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
      perc_paired <- -33
  }

  print("exited the dicer/fold segment.")
  if(exists("minus_null_res")){
    print("minus_null_res exists.")
    minus_res <- minus_null_res
  } else {
    print("minus_null_res does not exist.")
      minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
      minus_overhangz <- mean(minus_overhangs$zscore[1:4])
      minus_res <- list(minusMFE = MFE, minus_hp_overhangz = minus_hp_overhangz, minus_hp_phasedz = minus_hp_phasedz, phased_tbl.dist = minus_hp_phased_tbl$phased_dist,
                 phased_tbl.zscore = minus_hp_phased_tbl$phased_z, dicer_tbl.shift = minus_overhangs$shift, dicer_tbl.zscore = minus_overhangs$zscore, perc_paired = perc_paired)
  }
################################################################ make plots #####################################################################

  plus_overhang_out <- data.frame(t(plus_res$dicer_tbl.zscore))
  plus_overhang_out <- plus_overhang_out %>% dplyr::mutate(locus = paste0(chrom_name, "_", reg_start, "_", reg_stop))
  colnames(plus_overhang_out) <- c(plus_res$dicer_tbl.shift, 'locus')

  #if calc_overhang there are 10 columns
  plus_overhang_out <- plus_overhang_out[, c(10, 1:9)]

  #if calc_expand_overhang there are 18 columns
  #plus_overhang_out <- plus_overhang_out[, c(18, 1:17)]

  suppressWarnings(
     if(!file.exists("plus_hp_dicerz.txt")){
        write.table(plus_overhang_out, file = paste0(wkdir, "plus_hp_dicerz.txt"), sep = "\t", quote = FALSE, append = T, col.names = T, na = "NA", row.names = F)
     } else {
        write.table(plus_overhang_out, file = paste0(wkdir, "plus_hp_dicerz.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
     }
  )

  minus_overhang_out <- data.frame(t(minus_res$dicer_tbl.zscore))
  colnames(minus_overhang_out) <- minus_res$dicer_tbl.shift
  minus_overhang_out$locus <- paste0(chrom_name, "_", reg_start, "_", reg_stop)

  minus_overhang_out <- minus_overhang_out[, c(10, 1:9)]
  #minus_overhang_out <- minus_overhang_out[, c(18, 1:18)]

  suppressWarnings(
     if(!file.exists("minus_hp_dicerz.txt")){
        write.table(minus_overhang_out, file = paste0(wkdir, "minus_hp_dicerz.txt"), sep = "\t", quote = FALSE, append = T, col.names = T, na = "NA", row.names = F)
     } else {
        write.table(minus_overhang_out, file = paste0(wkdir, "minus_hp_dicerz.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
     }
  )

  prefix <- paste0(chrom_name, "_", reg_start, "_", reg_stop)

  plus_phased_out <- t(c(prefix, t(plus_res$phased_tbl.zscore)))
  minus_phased_out <- t(c(prefix, t(minus_res$phased_tbl.zscore)))


  suppressWarnings(
     if(!file.exists("plus_hp_phasedz.txt")){
        write.table(plus_phased_out, file = paste0(wkdir, "plus_hp_phasedz.txt"), sep = "\t", quote = FALSE, append = FALSE, col.names = F, na = "NA", row.names = F)
     } else {
        write.table(plus_phased_out, file = paste0(wkdir, "plus_hp_phasedz.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
     }
  )

  suppressWarnings(
     if(!file.exists("minus_hp_phasedz")){
        write.table(minus_phased_out, file = paste0(wkdir, "minus_hp_phasedz.txt"), sep = "\t", quote = FALSE, append = FALSE, col.names = F, na = "NA", row.names = F)
     } else {
        write.table(minus_phased_out, file = paste0(wkdir, "minus_hp_phasedz.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
    }
  )


  if(plot_output == 'T'){
    plus_overhangs <- data.frame(shift = plus_res$dicer_tbl.shift, zscore = plus_res$dicer_tbl.zscore)
    plus_overhangs$zscore[is.na(plus_overhangs$zscore)] <- 0

    minus_overhangs <- data.frame(shift = minus_res$dicer_tbl.shift, zscore = minus_res$dicer_tbl.zscore)
    minus_overhangs$zscore[is.na(minus_overhangs$zscore)] <- 0

    plus_overhang_plot <- plot_overhangz(plus_overhangs, "+")
    minus_overhang_plot <- plot_overhangz(minus_overhangs, "-")

    data <- read_densityBySize(bam_obj, chrom_name, reg_start, reg_stop, bam_file, wkdir)

    density_plot <- plot_density(data, reg_start, reg_stop)
    arc_plot <- plot_helix("helix.txt")

    plus_phasedz <- plot_hp_phasedz(plus_hp_phased_tbl, "+")

    minus_phasedz <- plot_hp_phasedz(minus_hp_phased_tbl, "-")


    ## plot genome annotations (optional)
    if(annotate_region == "T"){
      gtf_plot <- plot_gtf(gtf_file, chrom_name, reg_start, reg_stop)
      left <- cowplot::plot_grid(arc_plot, gtf_plot, density_plot, rel_widths = c(1,1,1), ncol = 1, align = "vh", axis = "lrtb")

      # Draw combined plot
      right <- cowplot::plot_grid(plus_overhang_plot, minus_overhang_plot, plus_phasedz, minus_phasedz, ncol = 1, align = "vh", axis = "l", rel_widths = c(1,1,1), rel_heights = c(1, 1, 1, 1))
   } else {
      left <- cowplot::plot_grid(arc_plot, NULL, density_plot, rel_widths = c(0.9,1,1), rel_heights = c(1,0.1,1), ncol = 1, align = "vh", axis = "lrtb")
      # Draw combined plot
      right <- cowplot::plot_grid(plus_overhang_plot, minus_overhang_plot,plus_phasedz, minus_phasedz, ncol = 1, rel_heights = c(1,1,1,1), rel_widths = c(1,1,1,1), align = "vh", axis = "l")
    }

    final_plot <- cowplot::plot_grid(left, NULL, right, ncol = 3, rel_heights = c(1,0.1,1), rel_widths = c(1,0.1, 1))

    prefix <- paste0(wkdir, chrom_name, "-", reg_start, "_", reg_stop, "_", strand)

    if(out_type == "png" || out_type == "PNG"){
      grDevices::png(file = paste0(prefix, "_hairpin_fold.png"), height = 9, width = 9, units = "in", res = 300)
      print(final_plot)
      grDevices::dev.off()
    } else {
      grDevices::pdf(file = paste0(prefix, "_hairpin_fold.pdf"), height = 9, width = 9)
      print(final_plot)
      grDevices::dev.off()
    }
  }

 #for machine learning / run_all
 return(list(minus_res, plus_res))
}
