#' run_all function
#' run all modules of the package and store metrics for ML
#' @param chrom_name a string
#' @param reg_stop an integer
#' @param reg_start an integer
#' @param chromosome an integer representing the chromosome number
#' @param length an integer
#' @param bam_file a string
#' @param roi a string
#' @param genome_file a string
#' @param min_read_count an integer
#' @param si_pal a string
#' @param pi_pal a string
#' @param plot_output a string, "T" or "F"
#' @param path_to_RNAfold a string
#' @param annotate_bed a string, "T" or "F"
#' @param weight_reads a string, "T", or "F"
#' @param bed_file a string
#' @return results

#' @export



new_run_all <- function(chrom_name, reg_start, reg_stop, chromosome, length, bam_file, roi, genome_file, min_read_count, si_pal, pi_pal,
                        plot_output, path_to_RNAfold, annotate_bed, weight_reads, bed_file){

  width <- pos <- start <- end <- NULL

  local_ml <- data.table::data.table(locus = numeric(1), locus_length = numeric(1), unique_read_bias = numeric(1),strand_bias = numeric(1),
  perc_GC = numeric(1), ave_size = numeric(1), perc_first_nucT = numeric(1), perc_A10 = numeric(1),
  highest_si_col = numeric(1), num_si_dicer_reads = numeric(1),si_dicerz = numeric(1), si_phasedz = numeric(1),MFE = numeric(1), hp_dicerz = numeric(1),
  mirnaMFE = numeric(1), mirna_dicerz = numeric(1), pingpong_col = numeric(1), max_pi_count = numeric(1),
  max_piz_overlap = numeric(1), pi_phasedz = numeric(1), pi_phased26z = numeric(1), mi_perc_paired = numeric(1), hp_perc_paired = numeric(1), overlapz = numeric(1),
  shap_p = numeric(1), auc = numeric(1))

  local_ml$locus <- paste0(chrom_name, ":", reg_start, "-", reg_stop)

  local_ml$locus_length <- reg_stop - reg_start

  print(local_ml$locus_length)
  prefix <- paste0(chrom_name, "_", reg_start, "-", reg_stop)
  ###############################
  # process bam input files

  if(!dir.exists('run_all/')) dir.create('run_all/')
  all_dir <- 'run_all/'

  if(!file.exists(paste0(all_dir, 'run_all_logfile.txt'))) file.create(paste0(all_dir, 'run_all_logfile.txt'))
  logfile <- paste0(all_dir, "run_all_logfile.txt")

  bam_obj <- OpenBamFile(bam_file)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)
  chr_name <- names(bam_header[['targets']])
  chr_length <- unname(bam_header[['targets']])
  bam_header <- NULL
  chromosome <- which(chr_name == chrom_name)

  cat(file = logfile, paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start, " reg_stop: ", reg_stop, "\n"), append = TRUE)
  cat(file = logfile, "Filtering forward and reverse reads by length\n", append = TRUE)

  chromP <- getChrPlus(bam_obj, chrom_name, reg_start, reg_stop)
  chromM <- getChrMinus(bam_obj, chrom_name, reg_start, reg_stop)

  forward_dt <- data.table::setDT(makeBamDF(chromP)) %>%
    subset(width <= 32 & width >= 15) %>%
    dplyr::mutate(start = pos, end = pos + width - 1) %>%
    dplyr::select(-c(pos))
  chromP <- NULL

  reverse_dt <- data.table::setDT(makeBamDF(chromM)) %>%
    subset(width <= 32 & width >= 15) %>%
    dplyr::mutate(start = pos, end = pos + width - 1) %>%
    dplyr::select(-c(pos))
  chromM <- NULL

  sizes <- data.frame(width = c(forward_dt$width, reverse_dt$width))

  set.seed(1234)
  sample <- sizes[sample(1:nrow(sizes)),]
  sample <- head(sample, 5000)
  print(head(sample))

  if(length(sample) > 3 & !(length(unique(sample) == 1))){
    local_ml$shap_p <- unlist(unname(shapiro.test(sample)))[2]
  } else {
    local_ml$shap_p <- 2
  }

  all_data <- rbind(forward_dt, reverse_dt)


#d <- density.default(all_data$width)
print("about to do the histogram")
if(nrow(all_data) > 1){
  m <- mean(all_data$width)
  std <- sqrt(var(all_data$width))

  if(!std == 0){
    bin_width <- KernSmooth::dpih(all_data$width, scalest = "stdev")

    nbins <- seq(min(all_data$width) - bin_width,
           max(all_data$width) + bin_width,
           by = bin_width)

    hist(all_data$width, density=20, breaks = 5, prob=TRUE,
     xlab="x-variable", ylim=c(0, 2),
    main="normal curve over histogram")

    curve <- curve(dnorm(x, mean=m, sd=std),
      col="darkblue", lwd=2, add=TRUE, yaxt="n")

    local_ml$auc <- sum(diff(curve$x) * (head(curve$y,-1)+tail(curve$y,-1)))/2
  } else {
      local_ml$auc <- 0
  }

} else {
   local_ml$auc <- -1
}


  print("got past the histogram")


  if(nrow(forward_dt) == 0 && nrow(reverse_dt) == 0) return()

  total_read_count <- nrow(forward_dt) + nrow(reverse_dt)

  unique_read_count <- nrow(forward_dt %>% dplyr::distinct(start, end)) + nrow(reverse_dt %>% dplyr::distinct(start,end))

  local_ml$unique_read_bias <- unique_read_count/total_read_count

  ###############################
  # get extra metrics for ML

  perc_plus <- nrow(forward_dt)/(nrow(forward_dt) + nrow(reverse_dt))
  perc_minus <- nrow(reverse_dt)/(nrow(reverse_dt) + nrow(forward_dt))

  #combine perc minus and perc plus into "strand bias"

  if(perc_plus > perc_minus){
     local_ml$strand_bias <- perc_plus
  } else {
     local_ml$strand_bias <- perc_minus
  }

  local_ml$perc_GC <- get_GC_content(forward_dt, reverse_dt)
  read_dist <- get_read_dist(bam_obj, chrom_name, reg_start, reg_stop)

  ave_size <- highest_sizes(read_dist)

  local_ml$ave_size <- ave_size

  cat(file = logfile, "Creating size plots\n", append = TRUE)

  if(!dir.exists('run_all/size_plots/')) dir.create('run_all/size_plots/')

  if(plot_output == 'T'){
    size_dir <- 'run_all/size_plots/'
    size_plots <- plot_sizes(read_dist)
  }

  local_ml$perc_first_nucT <- first_nuc_T(forward_dt, reverse_dt)
  local_ml$perc_A10 <- get_nuc_10(forward_dt, reverse_dt)

  max_sizes <- NULL
  read_dist <- NULL
  ###############################
  # run siRNA function
  cat(file = logfile, "Begin siRNA function\n", append = TRUE)
  if(!dir.exists('run_all/siRNA_dir/')) dir.create('run_all/siRNA_dir/')

  si_dir <- 'run_all/siRNA_dir/'
  si_log <- file.create('si_logfile.txt')


  si_res <- run_siRNA_function(chrom_name, reg_start, reg_stop, length, min_read_count, genome_file, bam_file, si_log, si_dir, si_pal, plot_output, path_to_RNAfold,
                           annotate_bed, weight_reads, bed_file)

  max_si_heat <- get_max_si_heat(si_res)

  local_ml$highest_si_col <- max_si_heat$highest_si_col
  local_ml$si_dicerz <- si_res$si_dicer$Z_score[9]
  local_ml$num_si_dicer_reads <- si_res[[2]]$proper_count[5]/total_read_count
  local_ml$hp_perc_paired <- max(unlist(unname(si_res[[3]][[1]][3])), unlist(unname(si_res[[3]][[2]][3])))

  perc_paired_file <- paste0(si_dir, "perc_paired.txt")
  col_status <- ifelse(exists_not_empty(perc_paired_file), FALSE, TRUE)
  write.table(local_ml$hp_perc_paired, perc_paired_file, quote = FALSE, append = TRUE, col.names = col_status)


  ## add in get phased results here

  #### get hairpin-specific results
  plus_phasedz <- unlist(unname(si_res[[3]][[1]][6]))
  if(plus_phasedz[1] != -33){
    plus_mean <- mean(plus_phasedz[1:4])
  } else {
    plus_mean <- -33
  }

   minus_phasedz <- unlist(unname(si_res[[3]][[2]][6]))
  if(minus_phasedz[1] != -33){

    minus_mean <- mean(minus_phasedz[1:4])
  } else {
    minus_mean <- -33
  }

  local_ml$si_phasedz <- max(plus_mean, minus_mean)

  local_ml$MFE <- min(unlist(unname(si_res[[3]][[1]][1])), unlist(unname(si_res[[3]][[2]][1])))

  local_ml$hp_dicerz <- max(unlist(unname(si_res[[3]][[1]][2])), unlist(unname(si_res[[3]][[2]][2])))
  #local_ml$hp_phasedz <- max(unlist(unname(si_res[[3]][[1]][3])), unlist(unname(si_res[[3]][[2]][3])))

  hp_dicerz_file <- paste0(si_dir, "hp_dicerz.txt")
  col_status <- ifelse(exists_not_empty(hp_dicerz_file), FALSE, TRUE)
  #write.table(local_ml$hp_dicerz, hp_dicerz_file, quote = FALSE, append = TRUE, col.names = col_status)

  #hp_phasedz_file <- paste0(si_dir, "hp_phasedz.txt")
  #col_status <- ifelse(exists_not_empty(hp_phasedz_file), FALSE, TRUE)
  #write.table(local_ml$hp_phasedz, hp_phasedz_file, quote = FALSE, append = TRUE, col.names = col_status)

  #print(paste0('hp_dicerz: ', local_ml$hp_dicerz))
  #print(paste0('hp_phasedz: ', local_ml$hp_phasedz))
  si_res <- NULL
  max_si_heat <- NULL
  ###############################
  # run miRNA function
  cat(file = logfile, "Begin miRNA function\n", append = TRUE)
  if(!dir.exists('run_all/miRNA_dir/')) dir.create('run_all/miRNA_dir/')
  miRNA_dir <- 'run_all/miRNA_dir/'
  mi_log <- file.create(paste0(miRNA_dir, 'mi_logfile.txt'))
  mi_res <- run_miRNA_function(chrom_name, reg_start, reg_stop, chromosome, length, "+", min_read_count, genome_file, bam_file, mi_log, miRNA_dir, plot_output, path_to_RNAfold, weight_reads)

  #Look at first result
  mi_res <- mi_res[[1]]
  mirnaMFE_plus <- mi_res$mfe

  pp_plus <- mi_res$perc_paired
  mirna_dicerz_plus <- mi_res$overhangs$zscore[5]

  if(mi_res$overhangs$zscore[1] != "NaN"){
    plus_overlapz <- mean(mi_res$overhangs$Z_score[20:23])
  } else {
    plus_overlapz <- NA
  }

  mi_res <- run_miRNA_function(chrom_name, reg_start, reg_stop, chromosome, length, "-", min_read_count, genome_file, bam_file, mi_log, miRNA_dir, plot_output, path_to_RNAfold, weight_reads)

  mi_res <- mi_res[[1]]
  mirnaMFE_minus <- mi_res$mfe
  pp_minus <- mi_res$perc_paired
  mirna_dicerz_minus <- mi_res$overhangs$zscore[5]

  if(mi_res$overhangs$zscore[1] != "NaN"){
    minus_overlapz <- mean(mi_res$overhangs$Z_score[20:23])
  } else {
    minus_overlapz <- NA
  }

  if(mirna_dicerz_plus == "NaN"){
    mirna_dicerz_plus <- -33
  }
  if(mirna_dicerz_minus == "NaN"){
    mirna_dicerz_minus <- -33
  }

  local_ml$mi_perc_paired <- max(pp_plus, pp_minus)
  local_ml$mirnaMFE <- min(mirnaMFE_plus, mirnaMFE_minus)
  local_ml$mirna_dicerz <- max(mirna_dicerz_plus, mirna_dicerz_minus)


  if(is.na(minus_overlapz) & !is.na(plus_overlapz)){
    local_ml$overlapz <- plus_overlapz
  } else if(!is.na(minus_overlapz) & is.na(plus_overlapz)){
    local_ml$overlapz <- minus_overlapz
  }  else {
      local_ml$overlapz <- max(minus_overlapz, plus_overlapz)
  }
  ################################
  # run piRNA function
  cat(file = logfile, "Begin piRNA function\n", append = TRUE)
  if(!dir.exists('run_all/piRNA_dir/')) dir.create('run_all/piRNA_dir/')

  piRNA_dir <- 'run_all/piRNA_dir/'
  pi_log <- file.create(paste0(piRNA_dir, 'pi_logfile.txt'))
  pi_res <- run_piRNA_function(chrom_name, reg_start, reg_stop, length, bam_file, genome_file, pi_log, piRNA_dir, pi_pal, plot_output = "F")

  if(!is.na(pi_res[[1]][1])){
    if(sum(pi_res[[1]] != 0)){
    max_pi_heat <- get_max_pi_heat(pi_res)

    local_ml$pingpong_col <- max_pi_heat$highest_pi_col

    local_ml$max_pi_count <- max_pi_heat$highest_pi_count/total_read_count

    local_ml$max_piz_overlap <- get_max_zscore(unlist(pi_res[[2]]$Z_score), unlist(pi_res[[2]]$Overlap))[[1]]
    piz_overlap_file <- paste0(piRNA_dir, "max_piz_overlap.txt")
    col_status <- ifelse(exists_not_empty(piz_overlap_file), FALSE, TRUE)
    write.table(local_ml$max_piz_overlap, piz_overlap_file, quote = FALSE, append = TRUE, col.names = col_status)
    }
  } else {

    local_ml$highest_pi_col <- NA

    local_ml$max_pi_count <- NA

    local_ml$max_piz_overlap <- NA
  }

  pi_res <- NULL
  max_pi_heat <- NULL
  ###############################
  # run phased_piRNA function
  cat(file = logfile, "Begin phased_piRNA function\n", append = TRUE)

  if(!dir.exists('run_all/phased_dir/')) dir.create('run_all/phased_dir/')
  phased_dir <- 'run_all/phased_dir/'
  phased_log <- file.create(paste0(phased_dir, 'phased_logfile.txt'))
  phased_res <- phased_piRNA_function("+", chrom_name, reg_start, reg_stop, bam_file, phased_log, phased_dir, plot_output= "F")

  phasedz_plus <- phased_res[1]

  phasedz26_plus <- phased_res[2]

  phased_res <- phased_piRNA_function("-", chrom_name, reg_start, reg_stop, bam_file, phased_log, phased_dir, plot_output= "F")

  phasedz_minus <- phased_res[1]
  phasedz26_minus <- phased_res[2]

  local_ml$phasedz <- max(phasedz_plus, phasedz_minus)
  local_ml$phased26z <- max(phasedz26_plus, phasedz26_minus)

  pi_phasedz_file <- paste0(phased_dir, "pi_phasedz.txt")
  col_status <- ifelse(exists_not_empty(pi_phasedz_file), FALSE, TRUE)
  write.table(local_ml$phasedz, pi_phasedz_file, quote = FALSE, append = TRUE, col.names = col_status)
  ####################################################################
  # add results to table
  tbl_pref <- strsplit(roi, "[.]")[[1]][1]
  tmp <- unlist(strsplit(bam_file, "[/]"))
  input_pref <- tmp[length(tmp)]
  input_pref2 <- strsplit(input_pref, "[.]")[[1]][1]

  tbl_name <- paste0(tbl_pref, "_", input_pref2)
  df <- as.matrix(local_ml)
  print("writing to table")

  cat(file = logfile, "Writing results to table\n", append = TRUE)

  ml_file <- paste0(tbl_name, "_ml.txt")
  col_status <- ifelse(exists_not_empty(ml_file), FALSE, TRUE)
  print(paste0("col_status: ", col_status))
  utils::write.table(df, ml_file, sep = "\t", quote = FALSE, append = T, col.names = col_status, na = "NA", row.names = F)

  print("file has been written.")
}
