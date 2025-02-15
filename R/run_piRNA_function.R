#' Run the piRNA function
#' Finds overlapping ping-piRNAs, then performs single-strand phasing calculation
#' Outputs a plot, a table of values for a heatmap, and a table of z-scores.
#' @param chrom_name a string
#' @param reg_start a whole number
#' @param reg_stop a whole number
#' @param length an integer
#' @param bam_file a string
#' @param genome_file a string
#' @param logfile a string
#' @param wkdir a string
#' @param pal a string
#' @param plot_output a bool, TRUE or FALSE, default = TRUE
#' @param weight_reads a string, Determines whether read counts will be weighted and with which method. Valid options are "weight_by_prop", "locus_norm", a user-defined value, or "none". See MiSiPi documentation for descriptions of the weighting methods.
#' @param write_fastas a bool, Determines whether piRNA pairs will be written to fasta. Expected values are TRUE or FALSE
#' @param out_type The type of file to write the plots to. Options are "png" or "pdf". Default is PDF.
#' @return plots, heat results, and zdf
#' @export
#'
run_piRNA_function <- function(chrom_name, reg_start, reg_stop, length, bam_file, genome_file, logfile, wkdir, pal,
                               plot_output, weight_reads, write_fastas, out_type) {
  
  prefix <- .get_region_string(chrom_name, reg_start, reg_stop)
  width <- pos <- NULL
  bam_obj <- .open_bam(bam_file, logfile)

  # Get the reads from the BAM using Rsamtools
  print("Making chromP and chromM")
  cat(file = paste0(wkdir, logfile), "Making chromP and chromM\n", append = TRUE)
  chromP <- .get_chr(bam_obj, chrom_name, reg_start, reg_stop, strand = "plus")
  chromM <- .get_chr(bam_obj, chrom_name, reg_start, reg_stop, strand = "minus")

  print("Finished making chromP and chromM. Filtering forward and reverse reads by length.")

  ################################################################# ping pong piRNA ##############################################################
  cat(file = paste0(wkdir, logfile), paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start - 1, " reg_stop: ", reg_stop - 1, "\n"), append = TRUE)
  cat(file = paste0(wkdir, logfile), paste0("Filtering forward and reverse reads by length", "\n"), append = TRUE)


  # Make forward and reverse dataframes filtered for width
  forward_dt <- data.table::setDT(.make_si_BamDF(chromP)) %>%
    subset(width <= 32 & width >= 16) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1)

  reverse_dt <- data.table::setDT(.make_si_BamDF(chromM)) %>%
    subset(width <= 32 & width >= 16) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1)


  # for the read size dist plot
  cat(file = paste0(wkdir, logfile), paste0("Getting read size distribution.", "\n"), append = TRUE)

  read_dist <- .get_read_size_dist(forward_dt, reverse_dt)

  # Summarize data frames into unique reads with a column of duplicates
  forward_dt <- forward_dt %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n())

  reverse_dt <- reverse_dt %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n())
  
  reverse_dt <- reverse_dt %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n())

  chromP <- NULL
  chromM <- NULL

  #for piRNA output tables
  size_dist <- dplyr::bind_rows(forward_dt, reverse_dt) %>%
    dplyr::group_by(width) %>%
    dplyr::summarize(count = sum(count))

  .output_readsize_dist(size_dist, prefix, wkdir, strand = NULL, type = "piRNA")
  size_dist <- NULL

  # Expand the data frames back to full size and weight reads in some cases
  print("Getting weighted data.frames.")
  #include "T" argument to return read sequences
  
  # Get expanded-weighted reads
  
  # locus_norm is being given a different locus_read_count here in piRNA vs in siRNA & dual_strand_hairpin
  # I think it is incorrect here, but I need confirmation as to what it should be
  # In siRNA it is given the sum of the number of reads in the forward or reverse set
  # Here it is being given the reads from both the forward and reverse or double the reverse
  
  # if (weight_reads == "Locus_norm" | weight_reads == "locus_norm") {
  #   print("normalize read count to locus")
  #   forward_dt <- locus_norm(forward_dt, sum(forward_dt$count, reverse_dt$count))
  #   reverse_dt <- locus_norm(reverse_dt, sum(reverse_dt$count, reverse_dt$count))
  # }
  
  # For now, I am changing this to be similar to siRNA

  locus_length <- reg_stop - reg_start + 1
  
  forward_dt <- .weight_reads(forward_dt, weight_reads, locus_length, sum(forward_dt$count))
  reverse_dt <- .weight_reads(reverse_dt, weight_reads, locus_length, sum(reverse_dt$count))
  
  print("Completed getting weighted dataframes.")


  #if there are both forward and reverse results
  if (!nrow(forward_dt) == 0 && !nrow(reverse_dt) == 0) {
    # Resummarize the reads for more efficient processing
    f_summarized <- forward_dt %>%
      dplyr::group_by_all() %>%
      dplyr::count()
    
    r_summarized <- reverse_dt %>%
      dplyr::group_by_all() %>%
      dplyr::count()
    
    #get the overlapping read pairs
    overlaps <- .find_overlaps(r_summarized, f_summarized)
    
    print("Completed find_overlaps.")

    #1/9/24 now make_BamDF returns sequence too, so read sequences can be extracted from that
    # ignoring reads with same start/stop but internal mismatches from output fasta
    if (write_fastas == TRUE) {
      #proper_overlaps <- overlaps %>% dplyr::filter(r1_end - r2_start == 10)
      #overlaps <- NULL
      proper_overlaps <- overlaps %>% dplyr::mutate(overlap = dplyr::case_when(r1_start > r2_start ~ (r2_end - r1_start), r1_start < r2_start ~ (r1_end - r2_start),
                                                                                 (r1_start >= r2_start & r1_end <= r2_end) ~ (r1_end - r1_start + 1),
                                                                                 (r2_start >= r1_start & r2_end <= r1_end) ~ (r2_end - r2_start + 1))) %>%
                                                                                 dplyr::filter(overlap == 10)
      rreads <- data.frame()
      freads <- data.frame()
      tmp <- rbind(forward_dt, reverse_dt)

      #need to fix this and make more efficient
      for(i in 1:nrow(proper_overlaps)){
        tmp_r <- tmp[which(tmp$start == proper_overlaps$r1_start[i] & tmp$end == proper_overlaps$r1_end[i]), ] %>%
          dplyr::distinct(start, end, .keep_all = TRUE)
        tmp_f <- tmp[which(tmp$start == proper_overlaps$r2_start[i] & tmp$end == proper_overlaps$r2_end[i]), ] %>%
          dplyr::distinct(start, end, .keep_all = TRUE)


        rreads <- rbind(rreads, tmp_r)
        freads <- rbind(freads, tmp_f)
      }


      rreads <- rreads %>% dplyr::rename("r1_start" = "start", "r1_end" = "end", "r1_seq" = seq) %>% dplyr::select(-c(width, first,rname))
      freads <- freads %>% dplyr::rename("r2_start" = "start", "r2_end" = "end", "r2_seq" = seq) %>% dplyr::select(-c(width, first,rname))

      paired_seqs <- cbind(rreads, freads)

      rreads <- NULL
      freads <- NULL

      paired_seqs <- paired_seqs %>%
        dplyr::mutate(read1_seq = paste0(">",chrom_name, ":", paired_seqs$r1_start, "-", paired_seqs$r1_end, " ", paired_seqs$r1_seq), read2_seq = paste0(">",chrom_name, ":", paired_seqs$r2_start, "-", paired_seqs$r2_end, " " , paired_seqs$r2_seq))

      fastas <- paired_seqs %>% dplyr::select(c(read1_seq, read2_seq)) %>%
        dplyr::transmute(col1 = paste0(read1_seq, ",", read2_seq)) %>%
        tidyr::separate_rows(col1, sep = ",")

      fastas <- stringi::stri_split_regex(fastas$col1, " ")

      write.table(unlist(fastas), file = paste0(prefix, "_piRNA_pairs.fa"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

      paired_seqs <- NULL
      fastas <- NULL
    }

    cat(file = paste0(wkdir, logfile), paste0("Making counts table.", "\n"), append = TRUE)
    
    z_res <- make_count_table(forward_dt$start, forward_dt$end, forward_dt$width,
                              reverse_dt$start, reverse_dt$end, reverse_dt$width)

    cat(file = paste0(wkdir, logfile), paste0("Finding overlaps.", "\n"), append = TRUE)

    heat_results <- get_overlap_counts(forward_dt$start, forward_dt$end, forward_dt$width,
                                       reverse_dt$end, reverse_dt$start, reverse_dt$width,
                                       check_pi = TRUE)

    #new_heat_plot <- .plot_si_heat(new_heat, chrom_name, reg_start, reg_stop, wkdir, pal = pal)
    # calculate overlaps of all reads for output table
    overlaps <- overlaps %>%
      dplyr::mutate(overlap = dplyr::case_when(r1_start > r2_start ~ (r2_end - r1_start),
                                               r1_start < r2_start ~ (r1_end - r2_start),
                                               (r1_start >= r2_start & r1_end <= r2_end) ~ (r1_end - r1_start + 1),
                                               (r2_start >= r1_start & r2_end <= r1_end) ~ (r2_end - r2_start + 1)),
                    total_dupes = r1_dupes * r2_dupes)
    # call new_pi_overlaps


    overlap_out <- data.frame(t(z_res))
    colnames(overlap_out) <- overlap_out[1,]
    overlap_out <- overlap_out[-1,]

    suppressWarnings(
      if(!file.exists(paste0(wkdir, "piRNA_alloverlaps_counts.txt"))){
        utils::write.table(overlap_out, file = paste0(wkdir, "piRNA_alloverlaps_counts.txt"), sep = "\t", quote = FALSE, append = FALSE, col.names = TRUE, na = "NA", row.names = FALSE)
      } else {
        utils::write.table(overlap_out, file = paste0(wkdir, "piRNA_alloverlaps_counts.txt"), quote = FALSE, sep = "\t", col.names = FALSE, append = TRUE, na = "NA", row.names = FALSE)
    })

    #forward_dt <- NULL
    #reverse_dt <- NULL
    row.names(heat_results) <- c('16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32')

    colnames(heat_results) <- c('16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32')

    output <- t(c(prefix, as.vector(heat_results)))

    suppressWarnings(
    if(!file.exists(paste0(wkdir, "piRNA_heatmap.txt"))){
        utils::write.table(output, file = paste0(wkdir, "piRNA_heatmap.txt"), sep = "\t", quote = FALSE, append = FALSE, col.names = TRUE, na = "NA", row.names = F)
      } else {
        utils::write.table(output, file = paste0(wkdir, "piRNA_heatmap.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
      })
    output <- NULL

    # Put results into table
    z_df <- data.frame("Overlap" = z_res[ ,1], "Z_score" = .calc_zscore(z_res$count))
    #z_res <- NULL



  } else {
    z_df <- data.frame(Overlap = c(seq(4,30)), Z_score = c(rep(0, times = 27)))
    heat_results <- matrix(data = 0, nrow = 17, ncol = 17)
    row.names(heat_results) <- c('16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32')
    colnames(heat_results) <- c('16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32')
  }


  f_summarized <- NULL
  r_summarized <- NULL
  overlaps <- NULL
  
################################################################### phased piRNAs #########################################################################
  #prefix <- paste0(chrom_name, "_", reg_start, "-", reg_stop)

  ################## compute plus strand

  cat(file = paste0(wkdir, logfile), paste0("Running plus strand for phased piRNAs.", "\n"), append = TRUE)
  print("Calculating phasing on plus strand.")
  
  # Processing unistrand, so make copy of original read df to transform
  filter_r1_dt <- forward_dt %>%
    dplyr::filter(first == "T") %>%
    dplyr::mutate(end = start + (width - 1) + 59) %>%
    dplyr::select(-first) %>%
    dplyr::group_by_all() %>%
    dplyr::count()

  filter_r2_dt <- forward_dt %>%
    dplyr::filter(first == "T") %>%
    dplyr::select(-first) %>%
    dplyr::group_by_all() %>%
    dplyr::count()

  all_table <- data.table::data.table(phased_dist=seq(0,50),
                                      phased_num=rep(0, 51))

  if (!nrow(filter_r1_dt) == 0) {
    phased_plus_counts <- .calc_phasing(filter_r1_dt, filter_r2_dt, 59)
  } else {
    # set null results for machine learning if no reads
    phased_plus_counts <- data.table::data.table(phased_dist = c(seq(0,50)),
                                                 phased_num = c(rep(0, times = 51)),
                                                 phased_z = NA)
    phased_plus_counts <- data.table::setDT(dplyr::full_join(phased_plus_counts,
                                                             all_table,
                                                             by = "phased_dist", "phased_num")) %>%
      dplyr::select(-phased_num.y) %>%
      dplyr::rename('phased_num' = phased_num.x)
    #phased_plus_counts[is.na(phased_plus_counts)] <- NA
  }

  cat(file = paste0(wkdir, logfile), paste0("Calculating plus strand phasing.", "\n"), append = TRUE)

  # for looking at only reads >= 26nt
  over_26_dt <- forward_dt %>%
    subset(width >= 26) %>%
    dplyr::filter(first =="T") %>%
    dplyr::mutate(end = start + (width - 1) + 59) %>%
    dplyr::select(-first) %>%
    dplyr::group_by_all() %>%
    dplyr::count()
  
  if (!nrow(over_26_dt) == 0) {
    phased_26_plus_counts <- .calc_phasing(over_26_dt, over_26_dt, 59)
  } else {
    phased_26_plus_counts <- data.table::data.table(phased_dist =  c(seq(0,50)),
                                                    phased_num = c(rep(0, times = 51)),
                                                    phased_z = NA)
    phased_26_plus_counts <- data.table::setDT(dplyr::full_join(phased_26_plus_counts,
                                                                all_table,
                                                                by = "phased_dist", "phased_num")) %>%
      dplyr::select(-phased_num.y) %>%
      dplyr::rename('phased_num' = phased_num.x)
  }

  phased_26_plus_counts <- phased_26_plus_counts %>%
    dplyr::rename(phased26_dist = phased_dist,
                  phased26_num = phased_num,
                  phased26_z = phased_z)

  #combine the results tables
  plus_df <- cbind(phased_plus_counts, phased_26_plus_counts)

  phased_plus_output <- phased_plus_counts %>%
    dplyr::select(phased_z)

  #transform table to one line for writing output
  phased_plus_output <- t(c(prefix, t(phased_plus_output)))

  phased26_plus_output <- phased_26_plus_counts %>%
    dplyr::select(phased26_z)
  phased26_plus_output <- t(c(prefix, t(phased26_plus_output)))


  suppressWarnings(
    if (!file.exists(paste0(wkdir,"phased_plus_piRNA_zscores.txt"))){
      write.table(phased_plus_output,
                  file = paste0(wkdir, "phased_plus_piRNA_zscores.txt"),
                  sep = "\t",
                  quote = FALSE,
                  append = FALSE,
                  col.names = T,
                  na = "NA",
                  row.names = F)
    } else {
      write.table(phased_plus_output,
                  file = paste0(wkdir, "phased_plus_piRNA_zscores.txt"),
                  quote = FALSE,
                  sep = "\t",
                  col.names = F,
                  append = TRUE,
                  na = "NA",
                  row.names = F)
    }
  )

  suppressWarnings(
    if (!file.exists(paste0(wkdir,"phased26_plus_piRNA_zscores.txt"))){
      write.table(phased26_plus_output,
                  file = paste0(wkdir, "phased26_plus_piRNA_zscores.txt"),
                  sep = "\t",
                  quote = FALSE,
                  append = FALSE,
                  col.names = T,
                  na = "NA",
                  row.names = F)
    } else {
      write.table(phased26_plus_output,
                  file = paste0(wkdir, "phased26_plus_piRNA_zscores.txt"),
                  quote = FALSE,
                  sep = "\t",
                  col.names = F,
                  append = TRUE,
                  na = "NA",
                  row.names = F)
    }
  )

  ################ run minus strand
  cat(file = paste0(wkdir, logfile), paste0("Running minus strand phasing.", "\n"), append = TRUE)
  print("Running minus strand phasing.")
  
  filter_r1_dt <- reverse_dt %>%
    dplyr::filter(first == "T") %>%
    dplyr::mutate(end = start + (width - 1) + 59) %>%
    dplyr::select(-first) %>%
    dplyr::group_by_all() %>%
    dplyr::count()
  
  filter_r2_dt <- reverse_dt %>%
    dplyr::filter(first == "T") %>%
    dplyr::select(-first) %>%
    dplyr::group_by_all() %>%
    dplyr::count()

  all_table <- data.table::data.table(phased_dist=seq(0,50),
                                      phased_num=rep(0, 51))

  cat(file = paste0(wkdir, logfile), paste0("Calculating minus strand phasing.", "\n"), append = TRUE)
  if (!nrow(filter_r1_dt) == 0) {
    phased_minus_counts <- .calc_phasing(filter_r1_dt, filter_r2_dt, 59)
  } else {
    phased_minus_counts <- data.table::data.table(phased_dist = c(seq(0,50)),
                                                  phased_num = c(rep(0, times = 51)),
                                                  phased_z = NA)
    phased_minus_counts <- data.table::setDT(dplyr::full_join(phased_minus_counts,
                                                              all_table,
                                                              by = "phased_dist", "phased_num")) %>%
      dplyr::select(-phased_num.y) %>%
      dplyr::rename('phased_num' = phased_num.x)
  }

  # for looking at only reads >= 26nt
  over_26_dt <- reverse_dt %>%
    subset(width >= 26) %>%
    dplyr::filter(first == "T") %>%
    dplyr::mutate(end = start + (width - 1) + 59) %>%
    dplyr::select(-first) %>%
    dplyr::group_by_all() %>%
    dplyr::count()

  if (!nrow(over_26_dt) == 0) {
    phased_26_minus_counts <- .calc_phasing(over_26_dt, over_26_dt, 59)
  } else {
    phased_26_minus_counts <- data.table::data.table(phased_dist = c(seq(0,50)),
                                                     phased_num = c(rep(0, times = 51)),
                                                     phased_z = NA)
    phased_26_minus_counts <- data.table::setDT(dplyr::full_join(phased_26_minus_counts,
                                                                 all_table,
                                                                 by = "phased_dist", "phased_num")) %>%
      dplyr::select(-phased_num.y) %>%
      dplyr::rename('phased_num' = phased_num.x)
  }

  #make the results data table

  phased_26_minus_counts <- phased_26_minus_counts %>%
    dplyr::rename(phased26_dist = phased_dist,
                  phased26_num = phased_num,
                  phased26_z = phased_z)

  minus_df <- cbind(phased_minus_counts, phased_26_minus_counts)

  phased_minus_output <- phased_minus_counts %>%
    dplyr::select(phased_z)
  phased_minus_output <- t(c(prefix, t(phased_minus_output)))

  phased26_minus_output <- phased_26_minus_counts %>%
    dplyr::select(phased26_z)
  phased26_minus_output <- t(c(prefix, t(phased26_minus_output)))

  all_phased <- data.frame(dist = phased_minus_counts$phased_dist)
  all_phased <- all_phased %>%
    dplyr::mutate(count = phased_plus_counts$phased_num + phased_minus_counts$phased_num)
  all_phased$zscore <- .calc_zscore(all_phased$count)
  tbl <- all_phased %>%
    dplyr::select(zscore)
  tbl <- as.data.frame(t(tbl))
  tbl$locus <- prefix
  colnames(tbl) <- c(all_phased$dist, "locus")
  tbl <- dplyr::select(tbl,52,1:51)

  suppressWarnings(
    if (!file.exists(paste0(wkdir, "phased_minus_piRNA_zscores.txt"))){
      write.table(phased_minus_output,
                  file = paste0(wkdir, "phased_minus_piRNA_zscores.txt"),
                  sep = "\t",
                  quote = FALSE,
                  append = FALSE,
                  col.names = F,
                  na = "NA",
                  row.names = F)
    } else {
      write.table(phased_minus_output,
                  file = paste0(wkdir, "phased_minus_piRNA_zscores.txt"),
                  quote = FALSE,
                  sep = "\t",
                  col.names = F,
                  append = TRUE,
                  na = "NA",
                  row.names = F)
    }
  )

  suppressWarnings(
    if (!file.exists(paste0(wkdir,"phased26_minus_piRNA_zscores.txt"))){
      write.table(phased26_minus_output,
                  file = paste0(wkdir, "phased26_minus_piRNA_zscores.txt"),
                  sep = "\t",
                  quote = FALSE,
                  append = FALSE,
                  col.names = F,
                  na = "NA",
                  row.names = F)
    } else {
      write.table(phased26_minus_output,
                  file = paste0(wkdir, "phased26_minus_piRNA_zscores.txt"),
                  quote = FALSE,
                  sep = "\t",
                  col.names = F,
                  append = TRUE,
                  na = "NA",
                  row.names = F)
    }
  )

  #all phased is the sum of zscores for both plus and minus strand
  suppressWarnings(
    if(!file.exists(paste0(wkdir,"all_phased_piRNA_zscores.txt"))){
      write.table(tbl,
                  file = paste0(wkdir, "all_phased_piRNA_zscores.txt"),
                  sep = "\t",
                  quote = FALSE,
                  append = FALSE,
                  col.names = T,
                  na = "NA",
                  row.names = F)
    } else {
      write.table(tbl,
                  file = paste0(wkdir, "all_phased_piRNA_zscores.txt"),
                  quote = FALSE,
                  sep = "\t",
                  col.names = F,
                  append = TRUE,
                  na = "NA",
                  row.names = F)
    }
  )

  print("Completed calculations. Making plots.")
  print(paste0("sum heat results: ", sum(heat_results)))

#################################################################################################
### make plots
#  if(!sum(heat_results) == 0 && plot_output == TRUE){

  if (plot_output == TRUE) {
    cat(file = paste0(wkdir, logfile), paste0("Generating plots.", "\n"), append = TRUE)
    ### ping pong plots
    read_dist <- .get_read_dist(bam_obj, chrom_name, reg_start, reg_stop)

    ## calculate read density by size
    data <- read_densityBySize(bam_obj, chrom_name, reg_start, reg_stop, bam_file, wkdir)

    z_df$Z_score[is.na(z_df$Z_score)] <- 0

    z <- .plot_overlapz(z_df)
    dist_plot <- plot_sizes(read_dist)

    if ((reg_stop - reg_start) > 7000) {
      density_plot <- .plot_large_density(data, reg_start, reg_stop)
    } else {
      density_plot <- .plot_density(data, reg_start, reg_stop)
    }
    data <- NULL
    dist_plot <- plot_sizes(read_dist)

    plus_df$phased_z[is.na(plus_df$phased_z)] <- 0
    plus_df$phased26_z[is.na(plus_df$phased26_z)] <- 0
    minus_df$phased_z[is.na(minus_df$phased_z)] <- 0
    minus_df$phased26_z[is.na(minus_df$phased26_z)] <- 0

    plus_phased_plot <- .plot_phasedz(plus_df, "+")
    minus_phased_plot <- .plot_phasedz(minus_df, "-")


    if (sum(heat_results) > 0) {
      options(scipen = 999)
      heat_plot <- .plot_si_heat(heat_results, chrom_name, reg_start, reg_stop, wkdir, pal = pal)

      top_left <- cowplot::plot_grid(dist_plot, NULL, ggplotify::as.grob(heat_plot), ncol = 1, rel_widths = c(0.8,1,1),
                                     rel_heights = c(0.8,0.1,1), align = "vh", axis = "lrtb")
      top_right <- cowplot::plot_grid(z, NULL, plus_phased_plot, NULL, minus_phased_plot, ncol = 1, rel_widths = c(1,1,1,1,1),
                                      rel_heights = c(1,0.1,1,0.1,1), align = "vh", axis = "lrtb")

      top <- cowplot::plot_grid(top_left, NULL, top_right, ncol = 3, rel_widths = c(1,0.1,1))

      bottom <- cowplot::plot_grid(density_plot)
      ## phased plots
      #left null right null bottom
      all_plot <- cowplot::plot_grid(top, NULL, bottom, nrow = 3,  ncol = 1, rel_widths = c(0.9,0.9,0.9),
                                     rel_heights = c(1,0.1,0.8), align = "vh", axis = "lrtb")
      #grDevices::pdf(file = paste0(wkdir, chrom_name,"_", reg_start,"-", reg_stop, "_pi-zscore.pdf"), height = 10, width = 14)

    } else {
      top_left <- cowplot::plot_grid(dist_plot, ncol = 1, rel_widths = c(1),
                                   rel_heights = c(0.8,0.1,1), align = "vh", axis = "lrtb")
      top_right <- cowplot::plot_grid(z, NULL, plus_phased_plot, NULL, minus_phased_plot, ncol = 1, rel_widths = c(1,1,1,1,1),
                                   rel_heights = c(1,0.1,1,0.1,1), align = "vh", axis = "lrtb")

      top <- cowplot::plot_grid(top_left, NULL, top_right, ncol = 3, rel_widths = c(1,0.1,1))

      bottom <- cowplot::plot_grid(density_plot)
      ## phased plots
                                                                                                                     #left null right null bottom
      all_plot <- cowplot::plot_grid(top, NULL, bottom, nrow = 3,  ncol = 1, rel_widths = c(0.9,0.9,0.9),
                                   rel_heights = c(1,0.1,0.8), align = "vh", axis = "lrtb")
     #grDevices::pdf(file = paste0(wkdir, chrom_name,"_", reg_start,"-", reg_stop, "_pi-zscore.pdf"), height = 10, width = 14)

    }
    
    if (out_type == "png" || out_type == "PNG") {
      grDevices::png(file = file.path(wkdir, paste0(prefix, "_pi-zscore.png")), width = 10, height = 11, bg = "white", units = "in", res = 300)
      print(all_plot)
      grDevices::dev.off()
    } else {
      grDevices::pdf(file = file.path(wkdir, paste0(prefix, "_pi-zscore.pdf")), width = 10, height = 11)
      print(all_plot)
      grDevices::dev.off()
    }
  }


  if (!is.na(sum(phased_plus_counts$phased_z))) {
    # get average zscore for first 4 distances (1-4nt)
    ave_plus_z <- mean(phased_plus_counts$phased_z[1:4])
  } else {
    ave_plus_z <- -33
  }

  if (!is.na(sum(phased_26_plus_counts$phased26_z))) {
    ave_plus_26z <- mean(phased_26_plus_counts$phased26_z[1:4])
  } else {
    ave_plus_26z <- -33
  }

  if (!is.na(sum(phased_minus_counts$phased_z))) {
    ave_minus_z <- mean(phased_minus_counts$phased_z[1:4])
  } else {
    ave_minus_z <- -33
  }

  if (!is.na(sum(phased_26_minus_counts$phased26_z))) {
    ave_minus_26z <- mean(phased_26_minus_counts$phased26_z[1:4])
  } else {
    ave_minus_26z <- -33
  }

  #return(c(ave_z, ave_26z))
  #cat(file = paste0(wkdir, logfile), paste0("Returning results for ML table.", "\n"), append = TRUE)
  #results for ML table
  return(list(heat_results = heat_results, z_df = z_df, phased_plus_z = ave_plus_z, phased_26plus_z = ave_plus_26z,
              phased_minus_z = ave_minus_z, phased_26minus_z = ave_minus_26z))
}

