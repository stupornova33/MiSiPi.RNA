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
#' @param plot_output a bool, TRUE or FALSE
#' @param path_to_RNAfold a string
#' @param path_to_RNAplot a string
#' @param annotate_region a bool, TRUE or FALSE
#' @param weight_reads a bool, TRUE or FALSE
#' @param gtf_file a string
#' @param write_fastas a bool, TRUE or FALSE. Default is FALSE
#' @param out_type Specifies whether file types for plots are png or pdf. Default is pdf.
#' @return results

#' @export



new_run_all <- function(chrom_name, reg_start, reg_stop,
                        chromosome, length, bam_file,
                        roi, genome_file, min_read_count,
                        si_pal, pi_pal, plot_output,
                        path_to_RNAfold, path_to_RNAplot,
                        annotate_region,
                        weight_reads, gtf_file,
                        write_fastas, out_type) {


  #wrapped <- function(chrom_name, reg_start, reg_stop,
  #                    chromosome, length, bam_file,
  #                    roi, genome_file, min_read_count,
  #                    si_pal, pi_pal, plot_output,
  #                    path_to_RNAfold, path_to_RNAplot,
  #                    annotate_region,
  #                    weight_reads, gtf_file,
  #                    write_fastas, out_type) {

  width <- pos <- start <- end <- NULL

  # create empty data table for results
  local_ml <- data.table::data.table(locus = numeric(1),
                                     locus_length = numeric(1),
                                     log_shap_p = numeric(1),
                                     auc = numeric(1),
                                     strand_bias = numeric(1),
                                     perc_GC = numeric(1),
                                     ave_size = numeric(1),
                                     perc_first_nucT = numeric(1),
                                     perc_A10 = numeric(1),
                                     highest_si_col = numeric(1),
                                     si_dicerz = numeric(1),
                                     num_si_dicer_reads = numeric(1),
                                     hp_perc_paired = numeric(1),
                                     hp_phasedz = numeric(1),
                                     hp_mfe = numeric(1),
                                     hp_dicerz = numeric(1),
                                     mi_perc_paired = numeric(1),
                                     mirna_dicerz = numeric(1),
                                     mirna_mfe = numeric(1),
                                     mirna_overlapz = numeric(1),
                                     pingpong_col = numeric(1),
                                     max_pi_count = numeric(1),
                                     max_piz_overlap = numeric(1),
                                     pi_phasedz = numeric(1),
                                     pi_phased26z = numeric(1))

  print("setting locus")
  local_ml$locus <- paste0(chrom_name, ":", reg_start, "-", reg_stop)

  local_ml$locus_length <- reg_stop - reg_start + 1

  print(local_ml$locus_length)
  print("setting prefix")
  prefix <- paste0(chrom_name, "_", reg_start, "-", reg_stop)

####################################################################### process bam input files #############################################################################
  #

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

  print("filtering forward and reverse dts")
  forward_dt <- data.table::setDT(make_si_BamDF(chromP)) %>%
    subset(width <= 32 & width >= 18) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1) %>%
    dplyr::group_by_all() %>%
    # get the number of times a read occurs
    dplyr::summarize(count = dplyr::n()) %>%
    na.omit()

  reverse_dt <- data.table::setDT(make_si_BamDF(chromM)) %>%
    subset(width <= 32 & width >= 18) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1) %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n()) %>%
    na.omit()

  size_dist <- dplyr::bind_rows(forward_dt, reverse_dt) %>%
    dplyr::group_by(width) %>%
    dplyr::summarize(count = sum(count))

  print("output_readsize_dist")
  output_readsize_dist(size_dist, prefix, all_dir, strand = NULL, type = "all")

  chromP <- NULL
  chromM <- NULL
  size_dist <- NULL

############################################################################ get extra metrics for ML ###################################################################

  # calculate
  # shap_p
  # auc
  # strand_bias
  # perc_GC
  # ave_size
  # perc_first_nucT
  # perc_A10
  print('Getting extra metrics for ML')

  sizes <- data.frame(width = c(forward_dt$width, reverse_dt$width))

  set.seed(1234)
  sample <- sizes[sample(1:nrow(sizes)),]
  sample <- head(sample, 5000)
  print(head(sample))
  sizes <- NULL

  if((length(sample) > 3) && !(length(unique(sample)) == 1)){
    local_ml$log_shap_p <- log10(as.numeric(unlist(unname(shapiro.test(sample)))[2]))
  } else {
    local_ml$log_shap_p <- 2
  }

  total_read_count <- sum(forward_dt$count) + sum(reverse_dt$count)

  unique_read_count <- nrow(forward_dt) + nrow(reverse_dt)

  local_ml$unique_read_bias <- unique_read_count/total_read_count

  if(nrow(forward_dt) > 0){
    forward_dt <- no_weight(forward_dt, as.character(chrom_name))
  } else {
    forward_dt <- forward_dt %>%
      dplyr::select(-c(count))
  }

  if(nrow(reverse_dt) > 0){
    reverse_dt <- no_weight(reverse_dt, as.character(chrom_name))
  } else {
    reverse_dt <- reverse_dt %>% dplyr::select(-c(count))
  }

  print("about to do historgram")
  if (nrow(forward_dt) + nrow(reverse_dt) > 1) {
    all_widths <- c(forward_dt$width, reverse_dt$width)

    m <- mean(all_widths)
    std <- sd(all_widths)

    if (!(std == 0)) {
      bin_width <- KernSmooth::dpih(all_widths, scalest = "stdev")

      nbins <- seq(min(all_widths) - bin_width,
                   max(all_widths) + bin_width,
                   by = bin_width)

      hist(all_widths, density = 20, breaks = 5, prob = TRUE,
           xlab = "Read size", ylim = c(0, 2),
           main = "Normal curve over histogram")

      curve <- curve(dnorm(x, mean = m, sd = std),
                     col = "darkblue", lwd = 2, add = TRUE, yaxt = "n")

      local_ml$auc <- sum(diff(curve$x) * (head(curve$y, -1) + tail(curve$y, -1))) / 2
    } else {
      local_ml$auc <- 0
    }
  } else {
    local_ml$auc <- -1
  }

  all_widths <- NULL

  print("got past the histogram")

  # TO DO: make a null_res to return
  if(nrow(forward_dt) == 0 && nrow(reverse_dt) == 0) return()

  perc_plus <- nrow(forward_dt)/(nrow(forward_dt) + nrow(reverse_dt))
  perc_minus <- nrow(reverse_dt)/(nrow(reverse_dt) + nrow(forward_dt))

  #combine perc minus and perc plus into "strand bias"

  if(perc_plus > perc_minus){
     local_ml$strand_bias <- perc_plus
  } else {
     local_ml$strand_bias <- perc_minus
  }

  all_seqs <- c(forward_dt$seq, reverse_dt$seq)
  local_ml$perc_GC <- get_GC_content(all_seqs)

  read_dist <- get_read_size_dist(forward_dt, reverse_dt)

  ave_size <- highest_sizes(read_dist)

  local_ml$ave_size <- ave_size

  cat(file = logfile, "Creating size plots\n", append = TRUE)

  if (!dir.exists('run_all/size_plots/')) {
    dir.create('run_all/size_plots/')
  }

  if(plot_output == TRUE){
    size_dir <- 'run_all/size_plots/'
    size_plots <- plot_sizes(read_dist)
  }

  local_ml$perc_first_nucT <- first_nuc_T(forward_dt, reverse_dt)

  all_nuc_10 <- all_seqs %>%
    stringr::str_sub(10, 10)
  all_nuc_10_A <- sum(all_nuc_10 == "A")
  local_ml$perc_A10 <- all_nuc_10_A / length(all_nuc_10)

  all_seqs <- NULL
  all_nuc_10 <- NULL
  all_nuc_10_A <- NULL
  max_sizes <- NULL
  read_dist <- NULL












############################################################################ run siRNA function #######################################################################
  # calculate
  # highest_si_col
  # si_dicerz
  # num_si_dicer_reads
  # hp_perc_paired

  cat(file = logfile, "Begin siRNA function\n", append = TRUE)
  if(!dir.exists('run_all/siRNA_dir/')) {
    dir.create('run_all/siRNA_dir/')
  }

  si_dir <- 'run_all/siRNA_dir/'
  si_log <- 'si_logfile.txt'

  if (!file.exists(file.path(si_dir, si_log))) {
    file.create(file.path(si_dir, si_log))
  }

  si_res <- run_siRNA_function(chrom_name, reg_start, reg_stop,
                               length, min_read_count, genome_file,
                               bam_file, si_log, si_dir,
                               si_pal, plot_output, path_to_RNAfold,
                               annotate_region, weight_reads, gtf_file,
                               write_fastas, out_type)

  max_si_heat <- get_max_si_heat(si_res)

  local_ml$highest_si_col <- max_si_heat$highest_si_col
  si_dicerz <- si_res$si_dicer$Z_score[5]

  if(is.na(si_dicerz)){
    local_ml$si_dicerz <- -33
  } else {
    local_ml$si_dicerz <- si_dicerz
  }

  plus_perc_paired <- si_res[[3]][[2]]$perc_paired
  minus_perc_paired <- si_res[[3]][[2]]$perc_paired
  # changed 3/25 to be RPM
  local_ml$num_si_dicer_reads <- (si_res[[2]]$proper_count[5]*1000000)/total_read_count
  local_ml$hp_perc_paired <- max(plus_perc_paired, minus_perc_paired)

  perc_paired_file <- paste0(si_dir, "perc_paired.txt")
  col_status <- ifelse(exists_not_empty(perc_paired_file), FALSE, TRUE)
  write.table(local_ml$hp_perc_paired, perc_paired_file, quote = FALSE, append = TRUE, col.names = col_status)

######################################################################### get hairpin-specific results ###############################################################

  # hp_phasedz [ maximum value of plus_phasedz[1:4] and minus_phasedz[1:4] ]
  # MFE change to hp_mfe
  # hp_dicerz [ maximum value of plus_dicerz and minus_dicerz]
  print("getting hairpin-specific results")
  #plus_phasedz <- unlist(unname(si_res[[3]][[2]][6]))
  plus_phasedz <- si_res[[3]][[2]]$phased_tbl.phased_z
  #if(!plus_phasedz[1] == "NaN" && !plus_phasedz[1] == -33){
  if(!is.na(plus_phasedz[1] && !plus_phasedz[1] == -33)){
    plus_mean <- mean(plus_phasedz[1:4])
  } else {
    plus_mean <- -33
  }

   #minus_phasedz <- unlist(unname(si_res[[3]][[1]][6]))
  minus_phasedz <- si_res[[3]][[1]]$phased_tbl.phased_z
  #if(!minus_phasedz[1] == "NaN" && !minus_phasedz[1] == -33){
   if(!is.na(minus_phasedz[1] && !minus_phasedz[1] == -33)){
    minus_mean <- mean(minus_phasedz[1:4])
  } else {
    minus_mean <- -33
  }

  local_ml$hp_phasedz <- max(plus_mean, minus_mean)
  local_ml$hp_mfe <- min(unlist(unname(si_res[[3]][[1]][1])), unlist(unname(si_res[[3]][[2]][1])))

  plus_dicerz <- si_res[[3]][[2]]$plus_hp_overhangz
  minus_dicerz <- si_res[[3]][[1]]$minus_hp_overhangz

  if(is.na(plus_dicerz)){
    plus_dicerz <- -33
  }

  if(is.na(minus_dicerz)){
    minus_dicerz <- -33
  }


  local_ml$hp_dicerz <- max(plus_dicerz, minus_dicerz)
  #local_ml$hp_phasedz <- max(unlist(unname(si_res[[3]][[1]][3])), unlist(unname(si_res[[3]][[2]][3])))

  hp_dicerz_file <- paste0(si_dir, "hp_dicerz.txt")
  col_status <- ifelse(exists_not_empty(hp_dicerz_file), FALSE, TRUE)
  #write.table(local_ml$hp_dicerz, hp_dicerz_file, quote = FALSE, append = TRUE, col.names = col_status)

  si_res <- NULL
  max_si_heat <- NULL

  ############################################################################# run miRNA function ####################################################################
  # mi_perc_paired
  # mirna_dicerz
  # mirna_mfe
  # mirna_overlapz

  cat(file = logfile, "Begin miRNA function\n", append = TRUE)
  if(!dir.exists('run_all/miRNA_dir/')) {
    dir.create('run_all/miRNA_dir/')
  }

  miRNA_dir <- 'run_all/miRNA_dir/'
  mi_log <- 'mi_logfile.txt'

  if (!file.exists(file.path(miRNA_dir, mi_log))) {
    file.create(paste0(miRNA_dir, mi_log))
  }

  mi_res <- new_miRNA_function(chrom_name, reg_start, reg_stop,
                               chromosome, length, "+",
                               min_read_count, genome_file, bam_file,
                               mi_log, miRNA_dir,
                               plot_output,
                               path_to_RNAfold,
                               path_to_RNAplot,
                               weight_reads,
                               write_fastas,
                               out_type)

  #Look at first result
  #mi_res <- mi_res[[1]]
  mirnaMFE_plus <- mi_res$mfe

  pp_plus <- mi_res$perc_paired
  mirna_dicerz_plus <- mi_res$overhangs$zscore[5]

  if(mi_res$overhangs$zscore[1] != "NaN"){
    plus_overlapz <- mean(mi_res$overhangs$Z_score[17:19])
  } else {
    plus_overlapz <- NA
  }

  mi_res <- new_miRNA_function(chrom_name, reg_start, reg_stop,
                               chromosome, length, "-",
                               min_read_count, genome_file, bam_file,
                               mi_log, miRNA_dir,
                               plot_output,
                               path_to_RNAfold,
                               path_to_RNAplot,
                               weight_reads,
                               write_fastas,
                               out_type)

  #mi_res <- mi_res[[1]]
  mirnaMFE_minus <- mi_res$mfe
  pp_minus <- mi_res$perc_paired
  mirna_dicerz_minus <- mi_res$overhangs$zscore[5]

  if(mi_res$overhangs$zscore[1] != "NaN"){
    minus_overlapz <- mean(mi_res$overhangs$Z_score[17:19])
  } else {
    minus_overlapz <- NA
  }

  if(mirna_dicerz_plus == "NaN"){
    mirna_dicerz_plus <- -33
  }
  if(mirna_dicerz_minus == "NaN"){
    mirna_dicerz_minus <- -33
  }

  if(is.na(mirnaMFE_minus) && !is.na(mirnaMFE_plus)){
    local_ml$mirna_mfe <- mirnaMFE_plus
  } else if(is.na(mirnaMFE_plus) && !is.na(mirnaMFE_minus)){
      local_ml$mirna_mfe <- mirnaMFE_minus
  } else if(is.na(mirnaMFE_minus) && is.na(mirnaMFE_plus)){
      local_ml$mirna_mfe <- 0
  } else {
      local_ml$mirna_mfe <- min(mirnaMFE_plus, mirnaMFE_minus)
  }
  local_ml$mi_perc_paired <- max(pp_plus, pp_minus)
  local_ml$mirna_dicerz <- max(mirna_dicerz_plus, mirna_dicerz_minus)


  if(is.na(minus_overlapz) && !is.na(plus_overlapz)){
    local_ml$mirna_overlapz <- plus_overlapz
  } else if(!is.na(minus_overlapz) && is.na(plus_overlapz)){
      local_ml$mirna_overlapz <- minus_overlapz
  } else if(is.na(minus_overlapz) && is.na(plus_overlapz)){
      local_ml$mirna_overlapz <- -33
  }  else {
      local_ml$mirna_overlapz <- max(minus_overlapz, plus_overlapz)
  }

############################################################################# run piRNA function ####################################################################
  # calculates pingpong_col
  # max_pi_count
  # max_piz_overlap

  cat(file = logfile, "Begin piRNA function\n", append = TRUE)
  if(!dir.exists('run_all/piRNA_dir/')) dir.create('run_all/piRNA_dir/')

  piRNA_dir <- 'run_all/piRNA_dir/'
  pi_log <- "pi_logfile.txt"
  if (!file.exists(file.path(piRNA_dir, pi_log))) {
    file.create(paste0(piRNA_dir, pi_log))
  }

  pi_res <- run_piRNA_function(chrom_name, reg_start, reg_stop,
                               length, bam_file, genome_file,
                               pi_log, piRNA_dir, pi_pal,
                               plot_output = FALSE,
                               weight_reads,
                               write_fastas,
                               out_type)


  #if(!is.na(pi_res[[1]])){
  #  max_pi_heat <- get_max_pi_heat(pi_res)

  #  local_ml$pingpong_col <- max_pi_heat$highest_pi_col

  #  local_ml$max_pi_count <- max_pi_heat$highest_pi_count/total_read_count

  #  local_ml$max_piz_overlap <- get_max_zscore(unlist(pi_res[[2]]$Z_score), unlist(pi_res[[2]]$Overlap))[[1]]
  #  piz_overlap_file <- paste0(piRNA_dir, "max_piz_overlap.txt")
  #  col_status <- ifelse(exists_not_empty(piz_overlap_file), FALSE, TRUE)
  #  write.table(local_ml$max_piz_overlap, piz_overlap_file, quote = FALSE, append = TRUE, col.names = col_status)
  #} else {
  #    local_ml$pingpong_col <- -33
  #    local_ml$max_pi_count <- -33
  #    local_ml$max_piz_overlap <- -33
  #}
  #pi_res <- NULL
  #max_pi_heat <- NULL

  if(sum(pi_res[[1]]) != 0){
    max_pi_heat <- get_max_pi_heat(pi_res)
    local_ml$pingpong_col <- max_pi_heat$highest_pi_col
    # changed pi_count to CPM
    local_ml$max_pi_count <- ((max_pi_heat$highest_pi_count)*1000000)/total_read_count
    local_ml$max_piz_overlap <- get_max_zscore(unlist(pi_res$z_df$Z_score), unlist(pi_res$z_df$Overlap))[[1]]
    piz_overlap_file <- paste0(piRNA_dir, "max_piz_overlap.txt")
    col_status <- ifelse(exists_not_empty(piz_overlap_file), FALSE, TRUE)
    write.table(local_ml$max_piz_overlap, piz_overlap_file, quote = FALSE, append = TRUE, col.names = col_status)
  } else {
    local_ml$pingpong_col <- -33
    local_ml$max_pi_count <- -33
    local_ml$max_piz_overlap <- -33
  }


  max_pi_heat <- NULL

## extract phasing results
  phasedz_plus <- pi_res$phased_plus_z
  phasedz26_plus <- pi_res$phased_26plus_z
  phasedz_minus <- pi_res$phased_minus_z
  phasedz26_minus <- pi_res$phased_26minus_z


  if(is.na(phasedz_minus) && !is.na(phasedz_plus)){
    local_ml$pi_phasedz <- phasedz_plus
  } else if(!is.na(phasedz_minus) && is.na(phasedz_plus)){
    local_ml$pi_phasedz <- phasedz_minus
  } else if(is.na(phasedz_minus) && is.na(phasedz_plus)){
    local_ml$pi_phasedz <- -33
  }  else {
    local_ml$pi_phasedz <- max(phasedz_plus, phasedz_minus)
  }

  if(is.na(phasedz26_minus) && !is.na(phasedz26_plus)){
    local_ml$pi_phased26z <- phasedz26_plus
  } else if(!is.na(phasedz26_minus) && is.na(phasedz26_plus)){
    local_ml$pi_phased26z <- phasedz26_minus
  } else if(is.na(phasedz26_minus) && is.na(phasedz26_plus)){
    local_ml$pi_phased26z <- -33
  }  else {
    local_ml$pi_phased26z <- max(phasedz26_plus, phasedz26_minus)
  }

  pi_phasedz_file <- paste0(piRNA_dir, "pi_phasedz.txt")
  col_status <- ifelse(exists_not_empty(pi_phasedz_file), FALSE, TRUE)
  write.table(local_ml$pi_phasedz, pi_phasedz_file, quote = FALSE, append = TRUE, col.names = col_status)




  pi_res <- NULL

####################################################################### add results to table ########################################################################

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
  #}

  #tryCatch(
  #  wrapped(chrom_name, reg_start, reg_stop,
  #          chromosome, length, bam_file,
  #          roi, genome_file, min_read_count,
  #          si_pal, pi_pal, plot_output,
  #          path_to_RNAfold, path_to_RNAplot,
  #          annotate_region,
  #          weight_reads, gtf_file,
  #          write_fastas, out_type),
  #  error = function(e) {
  #    message(paste("Uh oh boss:", conditionMessage(e)))
  #  }
  #)


}
