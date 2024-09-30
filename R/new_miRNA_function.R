#' new_miRNA_function
#' processes reads according to miRNA_algorithm, new with rnaplot function
#' returns plots
#' @param chrom_name a string
#' @param reg_start a whole number
#' @param reg_stop a whole number
#' @param chromosome chrom names extracted from bam file
#' @param length length of chromosome of interest
#' @param strand a character passed in, "+" or "-"
#' @param min_read_count a whole number
#' @param genome_file a fasta file of chrom sequences
#' @param bam_file a BAM file
#' @param logfile a string
#' @param wkdir a string
#' @param plot_output a bool, default = TRUE
#' @param path_to_RNAfold a string
#' @param path_to_RNAplot a string
#' @param write_fastas a bool, TRUE or FALSE
#' @param weight_reads Determines whether read counts will be weighted and with which method. Valid options are "weight_by_prop", "locus_norm", a user-defined value, or "none". See MiSiPi documentation for descriptions of the weighting methods.
#' @param out_type The type of file to write the plots to. Options are "png" or "pdf". Default is PDF.
#' @return plots
#' @export

new_miRNA_function <- function(chrom_name, reg_start, reg_stop, chromosome, length,
                               strand, min_read_count, genome_file, bam_file, logfile, wkdir, plot_output, path_to_RNAfold, path_to_RNAplot,
                               write_fastas, weight_reads, out_type){
  print(paste0(chrom_name, "-",reg_start,"-", reg_stop))
  cat(file = paste0(wkdir,logfile), paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start, " reg_stop: ", reg_stop, "\n"), append = TRUE)
  pos <- count <- count.x <- count.y <- end <- r1_end <- r1_start <- dist <- r2_end <- r2_start <- lstop <- lstart <- r1_seq <- loop_seq <- r2_seq <- start <- whole_seq <- width <- NULL

  # do not run locus if length is > 300 - not a miRNA. Also avoids issue where user provides coordinates of miRNA cluster.
  if(reg_stop - reg_start > 300){
    cat(file = paste0(wkdir, logfile), "length of region is greater than 300. \n", append = TRUE)
    return(null_mi_res())
  }

  bam_obj <- OpenBamFile(bam_file, logfile)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)
  chr_name <- names(bam_header[['targets']])
  chr_length <- unname(bam_header[['targets']])

  bam_header <- NULL

  #for the read size distribution plot
  chrom_m <- getChrMinus(bam_obj, chrom_name, reg_start, reg_stop)
  chrom_p <- getChrPlus(bam_obj, chrom_name, reg_start, reg_stop)

  read_dist <- get_read_dist(bam_obj, chrom_name, reg_start, reg_stop)

  # Moved this code block up so that the getChr functions don't have to be called again
  if(strand == "-"){
    chrom <- chrom_m
  } else {
    chrom <- chrom_p
  }
  chrom_m <- NULL
  chrom_p <- NULL
  which <- GenomicRanges::GRanges(seqnames=chrom_name, IRanges::IRanges(reg_start, reg_stop))

  if(strand == "-"){
    bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = TRUE), what=c('rname', 'pos', 'qwidth'), which=which)
  } else if(strand == "+"){
    bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = FALSE), what=c('rname', 'pos', 'qwidth'), which=which)
  } else{
    bam_scan <- Rsamtools::ScanBamParam(what=c('rname', 'pos', 'strand','qwidth'), which=which)
  }


  ########################################################## main logic ################################################################
  ## make the read data tables

  # Processing one strand, so make copy of df to transform
  # 1/9/24 added seq back in for writing miRNA pairs (arms)

  chrom_df <- make_si_BamDF(chrom) %>%
    subset(width <= 32 & width >= 18) %>%
    dplyr::mutate(start = pos, end = pos + width - 1) %>%
    dplyr::select(-c(pos)) %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n())
  filter_r2_dt <- chrom_df %>% dplyr::mutate(rname = chrom_name)

  if(nrow(filter_r2_dt) == 0){
    return(null_mi_res())
  } else {
    if(weight_reads == "weight_by_prop"){
      r2_dt <- weight_by_prop(filter_r2_dt, as.character(chrom_name))
      r1_dt <- weight_by_prop(filter_r2_dt, as.character(chrom_name))
    } else if(weight_reads == "Locus_norm" | weight_reads == "locus_norm"){
      locus_length <- reg_stop - reg_start
      locus_read_count <- sum(filter_r2_dt$count)
      r2_dt <- locus_norm(filter_r2_dt, locus_read_count)
      r1_dt <- locus_norm(filter_r2_dt, locus_read_count)
    } else {
      r2_dt <- no_weight(filter_r2_dt, as.character(chrom_name))
      r1_dt <- no_weight(filter_r2_dt, as.character(chrom_name))
    }
    #transform ends of one set of reads
    r1_dt <- r1_dt %>% dplyr::mutate(end = end + 59)
    filter_r2_dt <- NULL
    if(nrow(r2_dt) == 0){
      return(null_mi_res())
    }
  }

  chrom <- NULL
  chrom_df <- NULL

  if(nrow(r1_dt) == 0 || nrow(r2_dt) == 0){
    cat(paste0(wkdir, logfile), "After filtering for width and strand, zero reads remain. Please check bam BAM file.\n", append = TRUE)
    return(null_mi_res())
  }

  pileup_start <- min(r1_dt$start)
  pileup_stop <- max(r2_dt$end)
  pileups <- get_read_pileups(pileup_start, pileup_stop, bam_scan, bam_file)
  #dt <- pileups %>% dplyr::group_by(pos) %>% dplyr::summarise(count = sum(count))


  empty_table <- data.frame(pos = c(seq(pileup_start, pileup_stop)), count = c(0))

  dt_table <- merge(empty_table, pileups, by = 'pos', all.x = TRUE) %>% dplyr::select(-c(count.x)) %>% dplyr::rename('count' = count.y)
  dt_table[is.na(dt_table)] = 0


  # returns overlaps
  # using find_overlaps instead of find_hp_overlaps
  # find_overlaps only requires one df passed in and doesn't transform end of reads to original
  # find_hp_overlaps requires two dfs and automatically transforms ends of reads


  r1_dt <- r1_dt[sample(1:nrow(r1_dt)),]
  r1_dt <- utils::head(r1_dt, 10000)

  r2_dt <- r2_dt[sample(1:nrow(r2_dt)),]
  r2_dt <- utils::head(r2_dt, 10000)

  #system.time(test_overlap <- new_find_overlaps(r2_dt))

  overlaps_tmp <- find_overlaps(r1_dt, r2_dt)

  overlaps <- overlaps_tmp %>%
    dplyr::mutate(r1_end = r1_end - 59) %>%
    dplyr::mutate(r1_width = r1_end - r1_start + 1)

  overlaps_tmp <- NULL


  #need to filter to remove self matches, get dist between end of R1 and start of R2
  size <- nrow(overlaps)
  overlaps <- get_nearby(overlaps$r1_start, overlaps$r1_end, overlaps$r2_start, overlaps$r2_end, 60, size)


  # only want read pairs where start position of read 2 is after the end of read 1
  overlaps <- overlaps %>%
    dplyr::filter(dist > 2) %>%
    #dplyr::mutate(end_r1 = start_r1 + widthx - 1, end_r2 = start_r2 + widthy - 1) %>%
    dplyr::select(-c(dist))

  # write out pairs here

  if(nrow(overlaps) == 0){
    cat(paste0(wkdir, logfile), "No overlapping reads found.\n", append = TRUE)
    return(null_mi_res())
  }


  read_pileups <- getPileupsMap(dt_table$pos, dt_table$count, overlaps$start_r1, overlaps$end_r1, overlaps$start_r2, overlaps$end_r2)# %>%

  overlaps <- NULL

  if(nrow(read_pileups) == 0){
    return(null_mi_res())
  }

  read_pileups <- read_pileups %>%
    dplyr::mutate(width_r1 = (r1_end - r1_start) + 1, width_r2 = (r2_end - r2_start) + 1)


  read_pileups <- read_pileups %>% dplyr::mutate("lstart" = r1_end + 1, "lstop" = r2_start - 1) %>%  #get rid of rows with loop length less than 3
    dplyr::filter((lstop - lstart + 1) > 2)

  #loop_coord <- loop_coord %>% dplyr::filter((lstop - lstart + 1) > 2)
  if(nrow(read_pileups) == 0){
    return(null_mi_res())
  }

  #remove results where loop sequence has greater than 5% of total read count
  total_count <- sum(pileups$count)


  #df <- loop_coord

  mygranges <- GenomicRanges::GRanges(
    seqnames = c(chrom_name),
    ranges = IRanges::IRanges(start=c(1), end=c(length)))
  #print(paste0("mygranges: ", mygranges))

  geno_seq <- Rsamtools::scanFa(genome_file, mygranges)
  geno_seq <- as.character(unlist(Biostrings::subseq(geno_seq, start = 1, end = length)))


  loop_seqs <- getFastas(geno_seq, read_pileups$lstart, read_pileups$lstop, nrow(read_pileups)) %>%
    dplyr::rename("lstart" = "start", "lstop" = "stop")

  r1_seqs <- getFastas(geno_seq, read_pileups$r1_start, read_pileups$r1_end, nrow(read_pileups)) %>%
    dplyr::rename("r1_start" = "start", "r1_end" = "stop")
  r2_seqs <- getFastas(geno_seq, read_pileups$r2_start, read_pileups$r2_end, nrow(read_pileups)) %>%
    dplyr::rename("r2_start" = "start", "r2_end" = "stop")

  read_pileups <- read_pileups %>% dplyr::mutate(w_start = r1_start, w_stop = r2_end)

  read_pileups$r1_seq <- r1_seqs$Seq
  read_pileups$r2_seq <- r2_seqs$Seq
  read_pileups$loop_seq <- loop_seqs$Seq

  read_pileups$whole_seq <- stringr::str_c(read_pileups$r1_seq, read_pileups$loop_seq, read_pileups$r2_seq)
  read_pileups <- read_pileups %>% dplyr::mutate(width = w_stop - w_start + 1) %>%
                     dplyr::distinct()

  #final_coord <- test


  # select unique combinations of r1_start/stop, r2_start/stop
  grouped <- read_pileups %>%
    dplyr::group_by(r1_start, r1_end, r2_start, r2_end) %>%
    dplyr::distinct() %>%
    dplyr::select(c(r1_start, r1_end, r1_count_avg, r2_start, r2_end, r2_count_avg)) %>%
    dplyr::rename("r1_alt_start" = "r1_start", "r1_alt_end" = "r1_end", "r2_alt_start" = "r2_start", "r2_alt_end" = "r2_end")

  grouped <- grouped %>% dplyr::mutate(Chrom = chrom_name, Reg_start = reg_start, Reg_stop = reg_stop) %>%
    dplyr::select(c(Chrom, Reg_start, Reg_stop, r1_alt_start, r1_alt_end, r1_count_avg, r2_alt_start, r2_alt_end))


  if(nrow(grouped > 1)){
    write.table(grouped, file = paste0(wkdir, "alt_miRNAs_coord.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    cat(file = paste0(wkdir,logfile), "Writing potential alternative miRNA start and stop coordinates to alt_miRNAs_coord.bed.", append = TRUE)
  }

  most_abundant_idx <- which((read_pileups$r1_count_avg + read_pileups$r2_count_avg) == max(read_pileups$r1_count_avg + read_pileups$r2_count_avg))
    most_abundant <- read_pileups[most_abundant_idx,]

  final <- most_abundant[1,]

  mx_idx <- which(c(final$r1_count_avg, final$r2_count_avg) == max(final$r1_count_avg, final$r2_count_avg))
  mx <- c(final$r1_count_avg, final$r2_count_avg)[mx_idx]


  if(length(mx_idx) < 2){
    #colors <- vector(mode = "character", length = 2)
    colors <- c(NA, NA)
    #set max value color to YELLOW
    colors[mx_idx] <- "GREEN"
    #set not max value to GREEN
    colors[is.na(colors)] <- "RED"
    #colors[which(colors != mx_idx)] <- "GREEN"
  } else {
    colors <- vector(mode = "character", length = 2)
    colors[1] <- "GREEN"
    colors[2] <- "GREEN"
  }

  # need to convert starts and stops of reads to 1-based for coloring structure

  diff <- final$r1_start - reg_start

  pos_df <- data.frame(r1_start = c(0), r1_end = c(0),
                       r2_start = c(0), r2_end = c(0))
  if(diff < 0){
    #if final start coord is less than zero, add diff (+1 because R is 1-based)
    pos_df$r1_start <- final$r1_start + abs(diff) - reg_start + 1
    pos_df$r1_end <- final$r1_end + abs(diff) - reg_start + 1
    pos_df$r2_start <- final$r2_start + abs(diff) - reg_start + 1
    pos_df$r2_end <- final$r2_end + abs(diff) - reg_start + 1
  } else {
    pos_df$r1_start <- final$r1_start - reg_start + 1
    pos_df$r1_end <- final$r1_end - reg_start + 1
    pos_df$r2_start <- final$r2_start - reg_start + 1
    pos_df$r2_end <- final$r2_end - reg_start + 1
  }


  #call rnaplot here
  final_seq <- final$whole_seq
  rna_plot(path_to_RNAfold, path_to_RNAplot, wkdir, pos_df, colors)
  #RNAfold wants U's not T's, so convert to U
  converted <- list(convertU(final_seq, 1))

  converted <- data.frame('V1' = unname(unlist(converted)))

  colnames(converted) <- paste0(">",chrom_name, "_", reg_start, "_", reg_stop)
  final$converted <- converted$V1

  write.table(converted, file = paste0(wkdir, "converted.fasta"), sep = "\n", append = FALSE, row.names = FALSE, quote = FALSE)

  geno_seq <- NULL


  print('making fold_list')
  ################################################################################################################
  #fold_short_rna folds a list of sequences whereas fold_long_rna only folds one
  fold_list <- mapply(fold_short_rna, final$w_start, final$w_stop, converted, path_to_RNAfold)

  #extracts the needed fields from the fold list, store in more compact list
  make_reduced_list <- function(x){
    start <- fold_list[[1]][[2]]
    stop <- fold_list[[1]][[3]]
    mfe <- fold_list[[1]][[1]]
    vienna <- fold_list[[1]][[5]]
    helix <- R4RNA::viennaToHelix(vienna)
    converted <- fold_list[[1]][[6]]
    extracted_df <- fold_list[[1]][[4]]
    return(list(start = start, stop = stop, mfe = mfe, vienna = vienna, helix = helix, converted = converted, extracted_df = extracted_df))
  }

  reduced_list <- lapply(1:length(fold_list), make_reduced_list)


  #make the plots for all the sequences in the "reduced_list"
 plot_rna_struct <- function(x){
    pos <- count <- count.x <- NULL
    prefix <- paste0(wkdir, chrom_name, "-", (reduced_list[[1]]$start - 1), "-", (reduced_list[[1]]$stop - 1))
    mfe <- reduced_list[[1]]$mfe
    print(prefix)

    # put the helix table into a data.frame
    helix_df <- data.frame(reduced_list[[1]]$helix)

    perc_paired <- (length(reduced_list[[1]]$helix$i)*2)/(reduced_list[[1]]$stop - reduced_list[[1]]$start)

    # transforms reads from one arm of hairpin to their paired position
    # makes a table of reads which are overlapping
    dicer_overlaps <- dicer_overlaps(r2_dt, reduced_list[[1]]$helix, chrom_name, reduced_list[[1]]$start)
    # summarize the counts by the # overlapping nucleotides
    z_res <- make_count_table(r1_dt$start, r1_dt$end, r1_dt$width, r2_dt$start, r2_dt$end, r2_dt$width)

    # make_count_table was originally written for piRNAs. Need to subtract 3 from each overlap size.
    z_res <- z_res %>% dplyr::mutate(overlap = overlap - 3)
    # create empty z_df
    z_df <- data.frame("Overlap" = z_res[ ,1], "Z_score" = calc_zscore(z_res$count))

    # calculate the zscores, if there are results
    if(is.na(dicer_overlaps[1,1]) | dicer_overlaps[1,1] == 0) {
      overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
      overhangs$zscore <- calc_zscore(overhangs$proper_count)
    } else {
      #if(write_fastas == TRUE) write_proper_overhangs(wkdir, prefix, overlaps, "_miRNA")
      print('making overhangs')
      overhangs <- data.frame(calc_overhangs(dicer_overlaps$r1_start, dicer_overlaps$r1_end,
                                             dicer_overlaps$r2_start, dicer_overlaps$r2_width))
      overhangs$zscore <- calc_zscore(overhangs$proper_count)
    }

    # transform data frame from table to row
    overhang_output <- data.frame(t(overhangs$zscore))
    colnames(overhang_output) <- overhangs$shift
    overhang_output$locus <- paste0(chrom_name, "_", reduced_list[[x]]$start, "_", reduced_list[[x]]$stop)
    overhang_output <- overhang_output[, c(10, 1:9)]

    dice_file <- paste0(wkdir, "miRNA_dicerz.txt")
    suppressWarnings(
      if(!file.exists(dice_file)){
        utils::write.table(overhang_output, file = dice_file, sep = "\t", quote = FALSE, append = T, col.names = T, na = "NA", row.names = F)
      } else {
        utils::write.table(overhang_output, file = dice_file, quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
      }
    )


    if(plot_output == TRUE){
      #make the plots
      dicer_sig <- plot_overhangz(overhangs, "+")
      print(overhangs)
      #make new pileups dt for structure

      # get the per-base coverage
      # returns a two column df with pos and coverage
      new_pileups <- get_read_pileups(reduced_list[[1]]$start, reduced_list[[1]]$stop, bam_scan, bam_file)  %>%
        dplyr::group_by(pos) %>% dplyr::summarise(count = sum(count))

      # make a table with the same positions but empty cov column for combining with pileups
      # necessary because the pileups table doesn't have all positions from the first nt to the last because
      # coverages of zero aren't reported
      # set these to zero below
      empty_table <- data.frame(pos = c(seq(reduced_list[[1]]$start, reduced_list[[1]]$stop)), count = c(0))



      density <- read_densityBySize(bam_obj, chrom_name, reg_start, reg_stop, bam_file, wkdir)
      density_plot <- plot_density(density, reg_start, reg_stop)

      dist_plot <- plot_sizes(read_dist)

      zplot <- plot_overlapz(z_df)

      if(reg_stop - reg_start <= 150){
        #if miRNA is long, use landscape orientation so the struct plot isn't squished.
        left_top <- cowplot::plot_grid(dist_plot, dicer_sig, ncol = 1, rel_widths = c(1,1), rel_heights = c(1,1), align = "vh", axis = "lrtb")
        right_top <- cowplot::plot_grid(NULL, density_plot,zplot, ncol = 1, rel_widths = c(1,1,1), rel_heights = c(0.4,1,1), align = "vh", axis = "lrtb")
        #bottom <- cowplot::plot_grid(struct_plot)
        #top <- cowplot::plot_grid(left_top, right_top, rel_heights = c(1,1), rel_widths = c(1,1), align = "vh", axis = "lrtb")
        #all_plot <- cowplot::plot_grid(top,NULL, bottom, rel_widths = c(1,1,1), ncol = 1, rel_heights = c(1,0.1,0.5), align = "vh", axis = "lrtb")

        all_plot <- cowplot::plot_grid(left_top, right_top, rel_heights = c(1,1), rel_widths = c(1,1), align = "vh", axis = "lrtb")
        if(out_type == "png" || out_type == "PNG"){
          grDevices::png(file = paste0(prefix,"_", strand, "_combined.png"), height = 11, width = 8, units = "in", res = 300)
          print(all_plot)
          grDevices::dev.off()
        } else {
          grDevices::pdf(file = paste0(prefix,"_", strand, "_combined.pdf"), height = 11, width = 8)
          print(all_plot)
        }

      } else {
        all_plot <- cowplot::plot_grid(dist_plot, dicer_sig, density_plot, zplot, ncol = 4, rel_widths = c(1,1,1.4,1), align = "vh", axis = "lrtb")
        #left_top <- cowplot::plot_grid(dist_plot, dicer_sig, ncol = 1, rel_widths = c(1,1), rel_heights = c(1,1), align = "vh", axis = "lrtb")
        #right_top <- cowplot::plot_grid(NULL, density_plot,zplot, ncol = 1, rel_widths = c(1,1,1), rel_heights = c(0.4,1,1), align = "vh", axis = "lrtb")
        #bottom <- cowplot::plot_grid(struct_plot)
        #top <- cowplot::plot_grid(left_top, right_top, rel_heights = c(1,1), rel_widths = c(1,1), align = "vh", axis = "lrtb")
        #all_plot <- cowplot::plot_grid(top,NULL, bottom, rel_widths = c(1,1,0.6), ncol = 1, rel_heights = c(1,0.1,0.8), align = "vh", axis = "lrtb")


        if(out_type == "png" || out_type == "png"){
          grDevices::png(file = paste0(prefix,"_", strand, "_combined.png"), height = 7, width = 15, units = "in", res = 300)
          print(all_plot)
          grDevices::dev.off()
        } else {
          grDevices::pdf(file = paste0(prefix,"_", strand, "_combined.pdf"), height = 7, width = 15)
          print(all_plot)
          grDevices::dev.off()
        }

      }

    }
    return(list("mfe" = mfe, "perc_paired" = perc_paired, "overhangs" = c(overhangs,z_df)))
  }

  results <- lapply(1:length(reduced_list), plot_rna_struct)

  return(results)
}
