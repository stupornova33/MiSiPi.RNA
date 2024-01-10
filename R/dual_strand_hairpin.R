#' function to run the hairpin algorithm
#' processes reads from bam object according to strand
#' plots the arc plot and read dist
#' @param chrom_name a string
#' @param reg_start a whole number
#' @param reg_stop a whole number
#' @param length a whole number
#' @param min_read_count a whole number
#' @param genome_file a fasta file of chrom sequences
#' @param bam_file a BAM file
#' @param logfile a string
#' @param wkdir a string
#' @param plot_output a string, 'T' or 'F', default = 'T
#' @param path_to_RNAfold a string
#' @param annotate_bed a string, "T" or "F"
#' @param weight_reads a string, "T" or "F"
#' @param bed_file a string
#' @return max_overhang

#' @export

dual_strand_hairpin <- function(chrom_name, reg_start, reg_stop, length,
                             min_read_count, genome_file, bam_file, logfile, wkdir, plot_output, path_to_RNAfold, annotate_bed,
                             weight_reads, bed_file){

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
    res <- null_hp_res()
    return(res)
  }
  bam_header <- NULL

  # Extract chromosome sequence from genome file
  mygranges <- GenomicRanges::GRanges(
    seqnames = c(chrom_name),
    ranges = IRanges::IRanges(start=c(1), end=c(length)))

  geno_seq <- Rsamtools::scanFa(genome_file, mygranges)
  geno_seq <- as.character(unlist(Biostrings::subseq(geno_seq, start = 1, end = length)))

  cat(file = paste0(wkdir, logfile), paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start, " reg_stop: ", reg_stop, "\n"), append = TRUE)
  cat(file = paste0(wkdir, logfile), "Filtering forward and reverse reads by length\n", append = TRUE)

  #define which for Rsamtools ScanBamParam
  which <- GenomicRanges::GRanges(seqnames=chrom_name, IRanges::IRanges(reg_start, reg_stop))

  ############################################################ compute plus strand ########################################################
  strand <- "+"
  bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = FALSE), what=c('rname', 'pos', 'qwidth'), which=which)
  chrom <- getChrPlus(bam_obj, chrom_name, reg_start, reg_stop)


  #filter the reads and calculate the end position
  filter_r2_dt <- data.table::setDT(makeBamDF(chrom)) %>%
    base::subset(width <= 32 & width >= 18) %>%
    dplyr::mutate(start = pos, end = pos + width - 1) %>%
    dplyr::select(-c(pos)) %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n())

  # We're operating on a single strand but need one data frame where the reads aren't transformed, one where they are
  # so set the other dt to be the same as the first
  dt <- filter_r2_dt

  # weight reads if argument supplied
  if(weight_reads == "Top"){
    r2_dt <- get_top_n_weighted(dt, chrom_name, NULL)
  } else if(weight_reads == "Locus_norm"){
    r2_dt <- locus_norm(dt, sum(dt$count), NULL)
  } else {
    r2_dt <- get_top_n(dt, chrom_name, NULL)
  }

  #transform end of reads in one df
  r1_dt <- r2_dt %>% dplyr::mutate(end = end + 30)

  # if no results, need to store a specific value for the run_all/machine learning
  # null_hp_res() creates a table of specific "no result" values for zscores and such
  if(nrow(r1_dt) < 3 || nrow(r2_dt) < 3){
    cat(file = paste0(wkdir, logfile), "After filtering for width and strand, zero reads remain. Please check input BAM file.\n", append = TRUE)
    res <- null_hp_res()
    return(res)
  }

  # calculate phasing signatures
  if(nrow(r1_dt) > 0 && nrow(r2_dt) > 0){
    plus_hp_phased_tbl <- calc_phasing(r1_dt, r2_dt, 30)
    plus_hp_phased_counts <- sum(plus_hp_phased_tbl$phased_num[1:4])
    plus_hp_phased_z <- mean(plus_hp_phased_tbl$phased_z[1:4])
  } else {
    # if read dfs are empty set results to null. Still need to create the empty tables for plots/ML
    cat(file = paste0(wkdir, logfile), "No overlapping reads detected on this strand.\n", append = TRUE)
    # creating an empty table with "null" values
    plus_hp_phased_tbl <- data.table::data.table(phased_dist = seq(0,50), phased_num = rep(0,51), phased_z = rep(0,51))
    plus_hp_phased_counts <- sum(plus_hp_phased_tbl$phased_num[1:4])
    # -33 is an arbitrary value
    plus_phased_hp_z <- -33
  }

  if(nrow(r2_dt) > 0){

    # don't want to fold the dna for each strand, since it is the same. So once it's been folded
    # set fold_bool to TRUE
    fold_bool <- 'TRUE'

    fold_list <- fold_the_rna(geno_seq, chrom_name, reg_start, reg_stop, path_to_RNAfold)
    MFE <- fold_list$MFE
    perc_paired <- (length(fold_list$helix$i)*2)/(reg_stop - reg_start)

    #transform reads and find dicer pairs
    all_overlaps <- dicer_overlaps(r2_dt, fold_list$helix, chrom_name, reg_start)

    # dicer_overlaps() returns zero values if there are no valid overlaps
    # so check to make sure the first values are not zero
    if(!is.na(all_overlaps[1,1]) && !(all_overlaps[1,1] == 0)){
      plus_overhangs <- data.frame(calc_expand_overhangs(all_overlaps$r1_start, all_overlaps$r1_end,
                                               all_overlaps$r2_start, all_overlaps$r2_width))
      plus_overhangs$zscore <- calc_zscore(plus_overhangs$proper_count)

    } else {

      # return arbitrary "null" value if there are no valid results for ML
      plus_hp_phased_counts <- 0
      plus_hp_phased_z <- -33

      plus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
      plus_overhangs$zscore <- calc_zscore(plus_overhangs$proper_count)
    }
  } else {
      plus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
      plus_overhangs$zscore <- calc_zscore(plus_overhangs$proper_count)
      plus_hp_phased_counts <- sum(plus_overhangs$proper_count[1:4])
      plus_hp_phased_z <- mean(plus_overhangs$zscore[1:4])

      # if there were no results, and the dna didn't get folded, set fold_bool to FALSE so it gets folded in the minus strand part
      fold_bool <- 'FALSE'
      perc_paired <- 0
  }


  # results for the ML table
  plus_overhangs$zscore <- calc_zscore(plus_overhangs$proper_count)
  plus_overhangz <- mean(plus_overhangs$zscore[1:4])
  plus_res <- c(plusMFE = MFE, plus_hp_overhangz = plus_overhangz, plus_hp_phasedz = plus_hp_phased_z, phased_tbl = plus_hp_phased_tbl,
                dicer_tbl = plus_overhangs, perc_paired= perc_paired)

  ############################################################# compute minus strand ############################################################
  # do the same thing for the minus strand

  strand <- "-"
  bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = TRUE), what=c('rname', 'pos', 'qwidth'), which=which)
  chrom <- getChrMinus(bam_obj, chrom_name, reg_start, reg_stop)

  filter_r2_dt <- data.table::setDT(makeBamDF(chrom)) %>%
    base::subset(width <= 32 & width >= 18) %>%
    dplyr::mutate(start = pos, end = pos + width - 1) %>%
    dplyr::select(-c(pos))  %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n())


  if(weight_reads == "Top"){
    r2_dt <- get_top_n_weighted(filter_r2_dt, chrom_name, NULL)
  } else if(weight_reads == "Locus_norm"){
    r2_dt <- locus_norm(filter_r2_dt, sum(filter_r2_dt$count), NULL)
  } else {
    r2_dt <- get_top_n(filter_r2_dt, chrom_name, NULL)
  }

  #transform ends of reads for phasing/finding overlaps
  r1_dt <- r2_dt %>% dplyr::mutate(end = end + 30)

  if(nrow(r1_dt) < 3 || nrow(r2_dt) < 3){
    cat(file = paste0(wkdir, logfile), "After filtering for width and strand, zero reads remain. Please check input BAM file.\n", append = TRUE)
    res <- null_hp_res()
    return(res)
  }

  #calculate phasing
  if(nrow(r1_dt) > 0 && nrow(r2_dt) > 0){
     minus_hp_phased_tbl <- calc_phasing(r1_dt, r2_dt, 30)
     minus_hp_phased_counts <- sum(minus_hp_phased_tbl$phased_num[1:4])
     minus_phased_hp_z <- mean(minus_hp_phased_tbl$phased_z[1:4])
  } else {
    cat(file = paste0(wkdir, logfile), "No overlapping reads detected on this strand.\n", append = TRUE)
    minus_hp_phased_tbl <- data.table::data.table(phased_dist = seq(0,50), phased_num = rep(0,51), phased_z = rep(0,51))
    minus_hp_phased_counts <- sum(minus_hp_phased_tbl$phased_num[1:4])
    minus_phased_hp_z <- -33
  }

  if(nrow(r2_dt) > 0){
    if(fold_bool == 'TRUE'){
       all_overlaps <- dicer_overlaps(r2_dt, fold_list$helix, chrom_name, reg_start)

       if(!is.na(all_overlaps[1,1]) && !(all_overlaps[1,1] == 0)){  #if there are overlaps calc overhangs
         minus_overhangs <- data.frame(calc_expand_overhangs(all_overlaps$r1_start, all_overlaps$r1_end,
                                                     all_overlaps$r2_start, all_overlaps$r2_width))
         minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)

       } else {
           minus_hp_phased_counts <- 0
           minus_phased_hp_z <- -33

           minus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
           minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
           minus_hp_phased_counts <- sum(minus_overhangs$proper_count[1:4])
           minus_hp_phased_z <- mean(minus_overhangs$zscore[1:4])
       }
    } else { #else if fold bool is false and r2_dt > 0
        fold_list <- fold_the_rna(geno_seq, chrom_name, reg_start, reg_stop, converted, path_to_RNAfold)
        MFE <- fold_list$MFE
        perc_paired <- (length(fold_list$helix$i)*2)/(reg_stop - reg_start)
        all_overlaps <- dicer_overlaps(r2_dt, fold_list$helix, chrom_name, reg_start)

      if(!is.na(all_overlaps[1,1]) && !(all_overlaps[1,1] == 0)){  #if there are overlaps calc overhangs
         minus_overhangs <- data.frame(calc_expand_overhangs(all_overlaps$r1_start, all_overlaps$r1_end,
                                                      all_overlaps$r2_start, all_overlaps$r2_width))
         minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
      } else {
         minus_hp_phased_counts <- 0
         minus_phased_hp_z <- -33

         minus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
         minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
      }
   }

  } else { #else if fold bool is false and no results in r2_dt
      minus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
      minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
      perc_paired <- 0
  }

  minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
  minus_overhangz <- mean(minus_overhangs$zscore[1:4])
  minus_res <- c(minusMFE = MFE, minus_hp_overhangz = minus_overhangz, minus_hp_phasedz = minus_phased_hp_z, phased_tbl = minus_hp_phased_tbl,
                 dicer_tbl = minus_overhangs, perc_paired = perc_paired)

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
  colnames(minus_overhang_out) <- minus_res$shift
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

  plus_phased_out <- t(c(prefix, t(plus_res$phased_tbl.phased_z)))
  minus_phased_out <- t(c(prefix, t(minus_res$phased_tbl.phased_z)))


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
    minus_overhangs <- data.frame(shift = minus_res$dicer_tbl.shift, zscore = minus_res$dicer_tbl.zscore)

    plus_overhang_plot <- plot_overhangz(plus_overhangs, "+")
    minus_overhang_plot <- plot_overhangz(minus_overhangs, "-")

    data <- read_densityBySize(bam_obj, chrom_name, reg_start, reg_stop, bam_file, wkdir)

    density_plot <- plot_density(data, reg_start, reg_stop)
    arc_plot <- plot_helix("helix.txt")

    plus_phasedz <- plot_hp_phasedz(plus_hp_phased_tbl, "+")
    minus_phasedz <- plot_hp_phasedz(minus_hp_phased_tbl, "-")

    ## plot bed annotations (optional)
    if(annotate_bed == "T"){
      bed_plot <- plot_bed(bed_file, chrom_name, reg_start, reg_stop)
      left <- cowplot::plot_grid(arc_plot, bed_plot, density_plot, rel_widths = c(1,1,1), ncol = 1, align = "vh", axis = "lrtb")

      # Draw combined plot
      right <- cowplot::plot_grid(plus_overhang_plot, minus_overhang_plot, plus_phasedz, minus_phasedz, ncol = 1, align = "vh", axis = "l", rel_widths = c(1,1,1), rel_heights = c(1, 1, 1, 1))
   } else {
      left <- cowplot::plot_grid(arc_plot, NULL, density_plot, rel_widths = c(1,0.3,1), ncol = 1, align = "vh", axis = "lrtb")
      # Draw combined plot
      right <- cowplot::plot_grid(overhang_plot, NULL,phased_zscore, ncol = 1, rel_heights = c(1,0.1,1), rel_widths = c(1,1,1), align = "vh", axis = "l")
    }

    final_plot <- cowplot::plot_grid(left, NULL, right, ncol = 3, rel_heights = c(1,0.1,1), rel_widths = c(1,0.1, 1))

    prefix <- paste0(wkdir, chrom_name, "-", reg_start, "_", reg_stop, "_", strand)

    grDevices::pdf(file = paste0(prefix, "_hairpin_fold.pdf"), height = 8.5, width = 8.5)
    print(final_plot)
    grDevices::dev.off()
  }

 #for machine learning / run_all
 return(list(minus_res, plus_res))
}
