#' function to run the hairpin algorithm
#' processes reads from bam object according to strand
#' plots the arc plot and read dist 
#' @param chrom_name a string
#' @param reg_start a whole number
#' @param reg_stop a whole number
#' @param length
#' @param min_read_count a whole number
#' @param genome_file a fasta file of chrom sequences
#' @param input_file a BAM file
#' @param logfile a string
#' @param dir a string
#' @param plot_output a string, 'T' or 'F', default = 'T
#' @param path_to_RNAfold a string 
#' @param annotate_bed a string, "T" or "F"
#' @param gff_file a string
#' @return max_overhang

#' @export

dual_strand_hairpin <- function(chrom_name, reg_start, reg_stop, length,
                             min_read_count, genome_file, input_file, logfile, dir, plot_output, path_to_RNAfold, annotate_bed, 
                             gff_file){
  
  neg_results <- function(){
    MFE <- 0
    plus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
    plus_overhangs$zscore <- calc_zscore(plus_overhangs$proper_count)
    plus_hp_phased_counts <- 0
    plus_phased_hp_z <- -33
    minus_overhangs <- plus_overhangs 
    minus_hp_phased_counts <- 0
    minus_phased_hp_z <- -33
    return(list(MFE, plus_overhangs, minus_overhangs, plus_hp_phased_counts, minus_hp_phased_counts, plus_phased_hp_z, minus_phased_hp_z))
  }
  bam_obj <- OpenBamFile(input_file, logfile)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)
  
  
  if(reg_stop - reg_start > 10000){
    res <- neg_results()
    return(res)
  }
  bam_header <- NULL
  
  mygranges <- GenomicRanges::GRanges(
    seqnames = c(chrom_name),
    ranges = IRanges::IRanges(start=c(1), end=c(length)))
  print('made mygranges')
  
  geno_seq <- Rsamtools::scanFa(genome_file, mygranges) 
  geno_seq <- as.character(unlist(Biostrings::subseq(geno_seq, start = 1, end = length)))
  
  cat(file = paste0(dir, logfile), paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start, " reg_stop: ", reg_stop, "\n"), append = TRUE)
  cat(file = paste0(dir, logfile), "Filtering forward and reverse reads by length\n", append = TRUE)
  
  which <- GenomicRanges::GRanges(seqnames=chrom_name, IRanges::IRanges(reg_start, reg_stop))
  
    
  
  ################################################ compute plus strand ########################################################
  strand <- "+"
  
  if(strand == "-"){
    bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = TRUE), what=c('rname', 'pos', 'qwidth'), which=which)
  } else {
    bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = FALSE), what=c('rname', 'pos', 'qwidth'), which=which)
  }
  
  chrom_p <- getChrPlus(bam_obj, chrom_name, reg_start, reg_stop)
    
  filter_r2_dt <- get_top_n(chrom_p, chrom_name, 10)
  filter_r1_dt <- filter_r2_dt %>% dplyr::mutate(end = end + 59)
  
  if(nrow(filter_r1_dt) < 3 || nrow(filter_r2_dt) < 3){
    
    cat(file = paste0(dir, logfile), "After filtering for width and strand, zero reads remain. Please check input BAM file.\n", append = TRUE)
    res <- neg_results()
    return(res)
  }
  
  
  if(!nrow(filter_r1_dt) == 0 && !nrow(filter_r2_dt) == 0){
    print('making overlaps')
    overlaps <- find_hp_overlaps(filter_r1_dt, filter_r2_dt) 
  } else {
    overlaps <- data.frame(one = c(NA), two = c(NA))
  }
  
  
  print('making phased_counts')
  if(!is.na(overlaps[1,1])){
    phased_counts <- overlaps %>%
      dplyr::group_by(dist) %>%
      dplyr::summarize(num= dplyr::n())
    
    table <- data.table::data.table(dist=seq(1,50), num=rep(0, 50))
    phased_counts <- data.table::setDT(dplyr::full_join(phased_counts, table, by = "dist", "num"))
    
    phased_counts[is.na(phased_counts)] <- 0
    phased_counts <- phased_counts %>% dplyr::select(-c(num.y))
    phased_counts$Zscore <- calc_zscore(phased_counts$num.x)
    plus_hp_phased_tbl <- phased_counts %>% dplyr::rename(phased_dist = dist, phased_num = num.x, phased_z = Zscore)
    plus_phased_hp_z <- mean(plus_hp_phased_tbl$phased_z[1:4])
    plus_hp_phased_counts <- sum(plus_hp_phased_tbl$phased_num[1:4])
  } else {
    
    cat(file = paste0(dir, logfile), "No overlapping reads detected on this strand.\n", append = TRUE)
    #return(NA)
    plus_hp_phased_tbl <- data.table::data.table(phased_dist = seq(1,50), phased_num = rep(0,50), phased_z = rep(0,50))
    plus_hp_phased_counts <- sum(plus_hp_phased_tbl$phased_num[1:4])
    
    plus_phased_hp_z <- -33  #??
  }
  
  print('making dna_vec and converted')
  dna_vec <- as.character(Biostrings::subseq(geno_seq, start = reg_start, end = reg_stop))
  geno_seq <- NULL
  
  converted <- convertU(dna_vec, 1)
  writeLines(converted, con = "converted.txt")
  dna_vec <- NULL
  
  print('making fold_list')
  fold_list <- mapply(fold_long_rna, chrom_name, reg_start, reg_stop, converted, path_to_RNAfold)
  fold_list <- t(fold_list)
  MFE <- unlist(unname(fold_list[,3]))
  vienna <- fold_list[,5]
  extracted_df <- fold_list[4][[1]]
  
  writeLines(as.character(vienna), con = "vienna.txt")
  
  prefix <- paste0(dir, chrom_name, "-", reg_start, "_", reg_stop, "_", strand)
  
  print("making helix plot")
  helix <- R4RNA::viennaToHelix(unlist(fold_list[,5]))
  R4RNA::writeHelix(helix, file = "helix.txt")
  
  
  print('making all overlaps')
  
  if(nrow(filter_r2_dt) > 0){   #if there's stuff in the filter_dt, run dicer_overlaps ################################
    all_overlaps <- dicer_overlaps(filter_r2_dt, helix, chrom_name, reg_start) 
    
    if(!is.na(all_overlaps[1,1]) && !(all_overlaps[1,1] == 0)){  #if there are overlaps calc overhangs
      print('making overhangs')
      plus_overhangs <- data.frame(calc_overhangs(all_overlaps$r1_start, all_overlaps$r1_end,
                                             all_overlaps$r2_start, all_overlaps$r2_width))
      plus_overhangs$zscore <- calc_zscore(plus_overhangs$proper_count)
      #dicer_count <- overhangs$proper_count[5]
    } else {
      plus_hp_phased_counts <- 0
      plus_phased_hp_z <- -33
      
      plus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
      plus_overhangs$zscore <- calc_zscore(plus_overhangs$proper_count)
    }
  } else {
    #dicer_count <- 0
    plus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
    plus_overhangs$zscore <- calc_zscore(plus_overhangs$proper_count)
  }
  
  
  plus_overhangs$zscore <- calc_zscore(plus_overhangs$proper_count)
  plus_res <- c(MFE, plus_overhangs, plus_hp_phased_counts, plus_phased_hp_z)
  
  
  ######################################## compute minus strand without refolding ##############################################
  strand <- "-"
 
  if(strand == "-"){
    bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = TRUE), what=c('rname', 'pos', 'qwidth'), which=which)
  } else {
    bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = FALSE), what=c('rname', 'pos', 'qwidth'), which=which)
  }
  chrom_m <- getChrMinus(bam_obj, chrom_name, reg_start, reg_stop)
 
  
  filter_r2_dt <- get_top_n(chrom_m, chrom_name, 10)
  filter_r1_dt <- filter_r2_dt %>% dplyr::mutate(end = end + 59)
  
  if(nrow(filter_r1_dt) < 3 || nrow(filter_r2_dt) < 3){
    
    cat(file = paste0(dir, logfile), "After filtering for width and strand, zero reads remain. Please check input BAM file.\n", append = TRUE)
    res <- neg_results()
    return(res)
  }
  
  
  if(!nrow(filter_r1_dt) == 0 && !nrow(filter_r2_dt) == 0){
    print('making overlaps')
    overlaps <- find_hp_overlaps(filter_r1_dt, filter_r2_dt) 
  } else {
    overlaps <- data.frame(one = c(NA), two = c(NA))
  }
  
  
  print('making phased_counts')
  if(!is.na(overlaps[1,1])){
    phased_counts <- overlaps %>%
      dplyr::group_by(dist) %>%
      dplyr::summarize(num= dplyr::n())
    phased_counts <- head(phased_counts, 50)
    table <- data.table::data.table(dist=seq(1,50), num=rep(0, 50))
    phased_counts <- data.table::setDT(dplyr::full_join(phased_counts, table, by = "dist", "num"))
    
    phased_counts[is.na(phased_counts)] <- 0
    phased_counts <- phased_counts %>% dplyr::select(-c(num.y))
    phased_counts$Zscore <- calc_zscore(phased_counts$num.x)
    minus_hp_phased_tbl <- phased_counts %>% dplyr::rename(phased_dist = dist, phased_num = num.x, phased_z = Zscore)
    minus_phased_hp_z <- mean(minus_hp_phased_tbl$phased_z[1:4])
    minus_hp_phased_counts <- sum(minus_hp_phased_tbl$phased_num[1:4])
  } else {
    
    cat(file = paste0(dir, logfile), "No overlapping reads detected on this strand.\n", append = TRUE)
    #return(NA)
    minus_hp_phased_tbl <- data.table::data.table(phased_dist = seq(1,50), phased_num = rep(0,50), phased_z = rep(0,50))
    minus_hp_phased_counts <- sum(minus_hp_phased_tbl$phased_num[1:4])
    
    minus_phased_hp_z <- -33  #??
  }

  if(nrow(filter_r2_dt) > 0){   #if there's stuff in the filter_dt, run dicer_overlaps ################################
    all_overlaps <- dicer_overlaps(filter_r2_dt, helix, chrom_name, reg_start) 
    
    if(!is.na(all_overlaps[1,1]) && !(all_overlaps[1,1] == 0)){  #if there are overlaps calc overhangs
      print('making overhangs')
      minus_overhangs <- data.frame(calc_overhangs(all_overlaps$r1_start, all_overlaps$r1_end,
                                                  all_overlaps$r2_start, all_overlaps$r2_width))
      minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
      #dicer_count <- overhangs$proper_count[5]
    } else {
      minus_hp_phased_counts <- 0
      minus_phased_hp_z <- -33
      
      minus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
      minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
    }
  } else {
    #dicer_count <- 0
    minus_overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
    minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
  }
  
  minus_overhangs$zscore <- calc_zscore(minus_overhangs$proper_count)
  minus_res <- c(MFE, minus_overhangs, minus_hp_phased_counts, minus_phased_hp_z)

  ################################################### make plots ###########################################################
  
  print('making plots')
  if(plot_output == 'T'){
    print('making overhang plot')
    overhang_plot <- plot_overhangz(plus_overhangs, minus_overhangs)
    print('making data')
    data <- read_densityBySize(bam_obj, chrom_name, reg_start, reg_stop, input_file, dir)
    print('making density plot')
    density_plot <- plot_density(data, reg_start, reg_stop)
    print('making arc plot')
    arc_plot <- plot_helix("helix.txt")
    
    print('making phased_zscore')
    phased_zscore <- plot_phasedz(plus_hp_phased_tbl, minus_hp_phased_tbl)
   
    ## plot bed annotations (optional)
    if(annotate_bed == "T"){
      print("annotate_bed == T")
       
      bed_plot <- plot_bed(gff_file, chrom_name, reg_start, reg_stop)
      print('making left')
      left <- cowplot::plot_grid(arc_plot, bed_plot, density_plot, rel_widths = c(1,1,1), ncol = 1, align = "vh", axis = "lrtb")
   
      # Draw combined plot
      print('making right')
      right <- cowplot::plot_grid(overhang_plot, phased_zscore, ncol = 1, align = "vh", axis = "l")
   } else {
      left <- cowplot::plot_grid(arc_plot, NULL, density_plot, rel_widths = c(1,0.3,1), ncol = 1, align = "vh", axis = "lrtb")
      # Draw combined plot
      right <- cowplot::plot_grid(overhang_plot, phased_zscore, ncol = 1, align = "vh", axis = "l")
    }
    
    final_plot <- cowplot::plot_grid(left, right, ncol = 2, align = "vh", axis = "l", rel_widths = c(1, 0.9))
    prefix <- paste0(dir, chrom_name, "-", reg_start, "_", reg_stop, "_", strand)
    pdf(file = paste0(prefix, "_hairpin_fold.pdf"), height = 7, width = 7.5)
    print(final_plot)
    dev.off()
  }
  print("We made it boss!")
 # overhangs$zscore <- calc_zscore(overhangs$proper_count)
  #return(list(MFE, overhangs, hp_phased_counts, phased_hp_z))
}
