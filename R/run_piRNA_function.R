#' run the piRNA function
#' processes forward and reverse reads according to piRNA algorithm
#' outputs plots, heat results, and zdf
#' @param chrom_name a string
#' @param reg_start a whole number
#' @param reg_stop a whole number
#' @param bam_file a string
#' @param logfile a string
#' @param dir a string
#' @param pal a string
#' @param plot_output a string, 'T' or 'F', default = 'T
#' @return plots, heat results, and zdf
#' @export
#'
run_piRNA_function <- function(chrom_name, reg_start, reg_stop, bam_file, logfile, dir, pal, plot_output){
  width <- pos <- NULL
  bam_obj <- OpenBamFile(bam_file, logfile)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)

  bam_header <- NULL

  chromP <- getChrPlus(bam_obj, chrom_name, reg_start, reg_stop)
  chromM <- getChrMinus(bam_obj, chrom_name, reg_start, reg_stop)


  ### ping pong piRNA
  cat(file = paste0(dir, logfile), paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start, " reg_stop: ", reg_stop, "\n"), append = TRUE)
  cat(file = paste0(dir, logfile), paste0("Filtering forward and reverse reads by length", "\n"), append = TRUE)

  forward_dt <- data.table::setDT(makeBamDF(chromP)) %>%
    subset(width <= 32 & width >= 15) %>%
    dplyr::mutate(start = pos, end = pos + width - 1) %>% dplyr::distinct() %>%
    dplyr::select(-c(pos))

  reverse_dt <- data.table::setDT(makeBamDF(chromM)) %>%
    subset(width <= 32 & width >= 15) %>%
    dplyr::mutate(start = pos, end = pos + width - 1) %>% dplyr::distinct() %>%
    dplyr::select(-c(pos))

  chromP <- chromM <- NULL

  #### if no forward reads are appropriate length delete df and print error message
  if (nrow(forward_dt) == 0 || nrow(reverse_dt) == 0){
    cat(file = paste0(dir, logfile), paste0("Zero forward reads of correct length detected", "\n"), append = TRUE)
    z_df <- NA
    heat_results <- NA
    return(list(heat_results, z_df))
  }




  if(!nrow(forward_dt == 0) && (!nrow(reverse_dt) == 0)){
    overlaps <- find_overlaps(forward_dt, reverse_dt)
    size <- nrow(overlaps)

    mygranges <- GenomicRanges::GRanges(
    seqnames = c(chrom_name),
    ranges = IRanges::IRanges(start=c(1), end=c(length)))

    geno_seq <- Rsamtools::scanFa(genome_file, mygranges)
    geno_seq <- as.character(unlist(Biostrings::subseq(geno_seq, start = 1, end = length)))
    proper_overlaps <- overlaps %>% dplyr::filter(r1_end - r2_start == 10)
    paired_seqs <- proper_overlaps %>%
      dplyr::mutate(r1_seq = paste0(">",chrom_name, ":", proper_overlaps$r1_start, "-", proper_overlaps$r1_end, " ", substr(geno_seq, proper_overlaps$r1_start, proper_overlaps$r1_end)), r2_seq = paste0(">",chrom_name, ":", proper_overlaps$r2_start, "-", proper_overlaps$r2_end, " " , substr(geno_seq, proper_overlaps$r2_start, proper_overlaps$r2_end)))



    paired_seqs <- paired_seqs %>% dplyr::transmute(col1 = paste0(r1_seq, ",", r2_seq)) %>% tidyr::separate_rows(col1, sep = ",")
    fastas <- stringi::stri_split_regex(paired_seqs$col1, " ")

    write.table(unlist(fastas), paste0(dir, "piRNA_pairs.fa"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }



  z_res <- make_count_table(forward_dt$start, forward_dt$end, forward_dt$width,
                            reverse_dt$start, reverse_dt$end, reverse_dt$width)


  heat_results <- get_pi_overlaps(forward_dt$start, forward_dt$end, forward_dt$width,
                                  reverse_dt$end, reverse_dt$start, reverse_dt$width)

  forward_dt <- reverse_dt <- NULL

  row.names(heat_results) <- c('15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32')
  colnames(heat_results) <- c('15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32')

  prefix <- paste0(chrom_name, "_", reg_start, "_", reg_stop)
  output <- t(c(prefix, as.vector(heat_results)))

  suppressWarnings(
    if(!file.exists("piRNA_heatmap.txt")){
      utils::write.table(output, file = paste0(dir, "piRNA_heatmap.txt"), sep = "\t", quote = FALSE, append = T, col.names = F, na = "NA", row.names = F)
    } else {
      utils::write.table(output, file = paste0(dir, "piRNA_heatmap.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
    }
  )
  output <- NULL

  z_df <- data.frame("Overlap" = z_res[ ,1], "Z_score" = calc_zscore(z_res$count))
  z_res <- NULL
  #write.table(z_df, file = paste0(chrom_name, "_", reg_start, "-", reg_stop, "_zscore.txt"), row.names = FALSE, quote = FALSE


  ### phased piRNAs
  prefix <- paste0(chrom_name, "_", reg_start, "-", reg_stop)

  #for the read size dist plot

  read_dist <- get_read_dist(bam_obj, chrom_name, reg_start, reg_stop)

  ################## run plus strand

  chrom <- data.table::setDT(getChrPlus(bam_obj, chrom_name, reg_start, reg_stop))

  filter_dt <- data.table::setDT(makeBamDF(chrom)) %>%
    base::subset(width <= 32 & width >= 18) %>%
    dplyr::mutate(start = pos, end = pos + (width -1)) %>%
    dplyr::select(-c(rname, pos)) %>%
    dplyr::distinct()

  filter_dt <- filter_dt %>% dplyr::mutate(rname = chrom_name)

  filter_r1_dt <- filter_dt %>%
    dplyr::filter(first == "T") %>%
    dplyr::mutate(end = start + (width - 1) + 59) %>%
    dplyr::select(-c(first))

  over_26_dt <- filter_dt %>% subset(width >= 26)
  over_26_dt <- over_26_dt %>% dplyr::filter(first == "T") %>%
    dplyr::mutate(end = start + (width - 1) + 59) %>%
    dplyr::select(-c(first))

  filter_r2_dt <- filter_dt %>% dplyr::filter(first == "T") %>%
    dplyr::select(-c(first))

  all_table <- data.table::data.table(phased_dist=seq(1,50), phased_num=rep(0, 50))


  if(!nrow(filter_r1_dt) == 0) {
    phased_plus_counts <- calc_phasing(filter_r1_dt, filter_r2_dt)
  } else {
    phased_plus_counts <- data.table::data.table(phased_dist = 1, phased_num = 0L, phased_z = 0)
    phased_plus_counts <- data.table::setDT(dplyr::full_join(phased_plus_counts, all_table, by = "phased_dist", "phased_num")) %>%
      dplyr::select(-c(phased_num.y)) %>% dplyr::rename('phased_num' = phased_num.x)
    phased_plus_counts[is.na(phased_plus_counts)] <- 0
  }

  if(!nrow(over_26_dt) == 0){
    phased_26_plus_counts <- calc_phasing(over_26_dt, over_26_dt)
  } else {
    phased_26_plus_counts <- data.table::data.table(phased_dist = 1, phased_num = 0L, phased_z = 0)
    phased_26_plus_counts <- data.table::setDT(dplyr::full_join(phased_26_plus_counts, all_table, by = "phased_dist", "phased_num")) %>%
      dplyr::select(-c(phased_num.y)) %>% dplyr::rename('phased_num' = phased_num.x)
    phased_26_plus_counts[is.na(phased_26_plus_counts)] <- 0
  }

  #make the results data table

  phased_26_plus_counts <- phased_26_plus_counts %>%
    dplyr::rename(phased26_dist = phased_dist, phased26_num = phased_num, phased26_z = phased_z)


  plus_df <- cbind(phased_plus_counts, phased_26_plus_counts)

  prefix <- paste0(chrom_name, "_", reg_start, "_", reg_stop)


  phased_plus_output <- phased_plus_counts %>% dplyr::select(c(phased_z))
  phased_plus_output <- t(c(prefix, t(phased_plus_output)))

  phased26_plus_output <- phased_26_plus_counts %>% dplyr::select(c(phased26_z))
  phased26_plus_output <- t(c(prefix, t(phased26_plus_output)))


  suppressWarnings(
    if(!file.exists("phased_plus_piRNA_zscores.txt")){
      utils::write.table(phased_plus_output, file = paste0(dir, "phased_plus_piRNA_zscores.txt"), sep = "\t", quote = FALSE, append = FALSE, col.names = F, na = "NA", row.names = F)
    } else {
      utils::write.table(phased_plus_output, file = paste0(dir, "phased_plus_piRNA_zscores.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
    }
  )

  suppressWarnings(
    if(!file.exists("phased26_plus_piRNA_zscores.txt")){
      utils::write.table(phased26_plus_output, file = paste0(dir, "phased26_plus_piRNA_zscores.txt"), sep = "\t", quote = FALSE, append = FALSE, col.names = F, na = "NA", row.names = F)
    } else {
      utils::write.table(phased26_plus_output, file = paste0(dir, "phased26_plus_piRNA_zscores.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
    }
  )

  ############# run minus strand
  chrom <- getChrMinus(bam_obj, chrom_name, reg_start, reg_stop)

  filter_dt <- data.table::setDT(makeBamDF(chrom)) %>%
    base::subset(width <= 32 & width >= 18) %>%
    dplyr::mutate(start = pos, end = pos + (width - 1)) %>%
    dplyr::select(-c(rname, pos)) %>%
    dplyr::distinct()

  filter_dt <- filter_dt %>% dplyr::mutate(rname = chrom_name)

  over_26_dt <- filter_dt %>% subset(width >= 26)

  over_26_dt <- over_26_dt %>% dplyr::filter(first == "T") %>%
    dplyr::mutate(end = start + (width - 1) + 59) %>%
    dplyr::select(-c(first))

  filter_r1_dt <- filter_dt %>% dplyr::filter(first == "T") %>%
    dplyr::mutate(end = start + (width - 1) + 59) %>%
    dplyr::select(-c(first))

  filter_r2_dt <- filter_dt %>% dplyr::filter(first == "T") %>%
    dplyr::select(-c(first))



  over_26_dt <- filter_dt %>% subset(width >= 26)
  over_26_dt <- over_26_dt %>% dplyr::filter(first == "T") %>%
    dplyr::mutate(end = start + (width - 1) + 59) %>%
    dplyr::select(-c(first))

  filter_r2_dt <- filter_dt %>% dplyr::filter(first == "T") %>%
    dplyr::select(-c(first))

  all_table <- data.table::data.table(phased_dist=seq(1,50), phased_num=rep(0, 50))


  if(!nrow(filter_r1_dt) == 0) {
    phased_minus_counts <- calc_phasing(filter_r1_dt, filter_r2_dt)
  } else {
    phased_minus_counts <- data.table::data.table(phased_dist = 1, phased_num = 0L, phased_z = 0)
    phased_minus_counts <- data.table::setDT(dplyr::full_join(phased_minus_counts, all_table, by = "phased_dist", "phased_num")) %>%
      dplyr::select(-c(phased_num.y)) %>% dplyr::rename('phased_num' = phased_num.x)
    phased_minus_counts[is.na(phased_minus_counts)] <- 0
  }

  if(!nrow(over_26_dt) == 0){
    phased_26_minus_counts <- calc_phasing(over_26_dt, over_26_dt)
  } else {
    phased_26_minus_counts <- data.table::data.table(phased_dist = 1, phased_num = 0L, phased_z = 0)
    phased_26_minus_counts <- data.table::setDT(dplyr::full_join(phased_26_minus_counts, all_table, by = "phased_dist", "phased_num")) %>%
      dplyr::select(-c(phased_num.y)) %>% dplyr::rename('phased_num' = phased_num.x)
    phased_26_minus_counts[is.na(phased_26_minus_counts)] <- 0
  }

  #make the results data table

  phased_26_minus_counts <- phased_26_minus_counts %>%
    dplyr::rename(phased26_dist = phased_dist, phased26_num = phased_num, phased26_z = phased_z)


  minus_df <- cbind(phased_minus_counts, phased_26_minus_counts)

  prefix <- paste0(chrom_name, "_", reg_start, "_", reg_stop)


  phased_minus_output <- phased_minus_counts %>% dplyr::select(c(phased_z))
  phased_minus_output <- t(c(prefix, t(phased_minus_output)))

  phased26_minus_output <- phased_26_minus_counts %>% dplyr::select(c(phased26_z))
  phased26_minus_output <- t(c(prefix, t(phased26_minus_output)))


  suppressWarnings(
    if(!file.exists("phased_minus_piRNA_zscores.txt")){
      utils::write.table(phased_minus_output, file = paste0(dir, "phased_minus_piRNA_zscores.txt"), sep = "\t", quote = FALSE, append = FALSE, col.names = F, na = "NA", row.names = F)
    } else {
      utils::write.table(phased_minus_output, file = paste0(dir, "phased_minus_piRNA_zscores.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
    }
  )

  suppressWarnings(
    if(!file.exists("phased26_minus_piRNA_zscores.txt")){
      utils::write.table(phased26_minus_output, file = paste0(dir, "phased26_minus_piRNA_zscores.txt"), sep = "\t", quote = FALSE, append = FALSE, col.names = F, na = "NA", row.names = F)
    } else {
      utils::write.table(phased26_minus_output, file = paste0(dir, "phased26_minus_piRNA_zscores.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
    }
  )


  #################################################################################################
  ## make plots

  if(sum(heat_results != 0) && plot_output == 'T'){
    ### ping pong plots
    read_dist <- get_read_dist(bam_obj, chrom_name, reg_start, reg_stop)

    ## calculate read density by size
    data <- read_densityBySize(bam_obj, chrom_name, reg_start, reg_stop, bam_file, dir)

    z <- plot_overlapz(z_df)
    dist_plot <- plot_sizes(read_dist)

    heat_plot <- plot_si_heat(heat_results, chrom_name, reg_start, reg_stop, dir, pal = pal)

    density_plot <- plot_density(data, reg_start, reg_stop)
    data <- NULL
    dist_plot <- plot_sizes(read_dist)
    plus_phased_plot <- plot_phasedz(plus_df, "+")
    minus_phased_plot <- plot_phasedz(minus_df, "-")


    top <- cowplot::plot_grid(dist_plot, NULL, density_plot, ncol = 3, rel_widths = c(1,0.1,1), align = "vh", axis = "lrtb")
    middle <- cowplot::plot_grid(ggplotify::as.grob(heat_plot), NULL, z, rel_widths = c(1,0.1,0.8), nrow = 1, ncol = 3, align = "vh", axis = "lrtb")
    bottom <- cowplot::plot_grid(plus_phased_plot, NULL, minus_phased_plot, rel_widths = c(1,0.1,1), nrow = 1, ncol = 3, align = "vh", axis = "lrtb")
    ## phased plots

    all_plot <- cowplot::plot_grid(top, NULL, middle, NULL,bottom, ncol = 1, rel_widths = c(1,1, 1, 1,0.8), rel_heights = c(1,0.1, 1.2, 0.1, 0.8))
    grDevices::pdf(file = paste0(dir, chrom_name,"_", reg_start,"-", reg_stop, "_pi-zscore.pdf"), height = 7, width = 7)
    print(all_plot)
    grDevices::dev.off()
  }


  if(!is.na(sum(phased_plus_counts$phased_z))){

    ave_plus_z <- mean(phased_plus_counts$phased_z[1:4])
  } else {
    ave_plus_z <- -30
  }


  if(!is.na(sum(phased_26_plus_counts$phased26_z))){
    ave_plus_26z <- mean(phased_26_plus_counts$phased26_z[1:4])
  } else {
    ave_plus_26z <- -30
  }


  if(!is.na(sum(phased_minus_counts$phased_z))){

    ave_minus_z <- mean(phased_minus_counts$phased_z[1:4])
  } else {
    ave_minus_z <- -30
  }


  if(!is.na(sum(phased_26_minus_counts$phased26_z))){
    ave_minus_26z <- mean(phased_26_minus_counts$phased26_z[1:4])
  } else {
    ave_minus_26z <- -30
  }


  #return(c(ave_z, ave_26z))

  #results for ML table
  return(list(heat_results, z_df, ave_plus_z, ave_plus_26z, ave_minus_z, ave_minus_26z))
}

