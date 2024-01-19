#' Function to process and set the arguments passed by the user for each function
#' @param roi The path to a BED file of loci of interest
#' @param bam_file The path to a BAM file
#' @param genome The path to a genome Fasta file
#' @param min_read_count An integer. Default is 1
#' @param plot_output Determines whether the program will output plots as PDFs. Expected input is "T" or "F".
#' @param path_to_RNAfold The full path to the RNAfold binary executable.
#' @param pi_pal The color palette to use for the piRNA heatmap plot. Valid options are "RdYlBl", "BlYel", "yelOrRed", "MagYel", and "Greens".
#' @param si_pal The color palette to use for the siRNA heatmap plot. Valid options are "RdYlBl", "BlYel", "yelOrRed", "MagYel", and "Greens".
#' @param annotate_region Determines whether the program will plot genomic features of interest found in the GTF annotation file. If "T", a GTF file must be provided as the "gtf_file" argument.
#' @param weight_reads Determines whether read counts will be weighted. Valid options are "Top", "locus_norm", or "none". See MiSiPi documentation for descriptions of the weighting methods.
#' @param gtf_file a string corresponding to the path of genome annotation in 9-column GTF format.
#' @param write_fastas A string, "T" or "F". Optional. If "T", read pairs from functions will be written to file.
#'
#' @return a list
#' @export

set_vars <- function(roi, bam_file, genome, min_read_count, plot_output, path_to_RNAfold, pi_pal, si_pal, annotate_region, weight_reads, gtf_file = NULL, write_fastas = NULL){

  bam_obj <- OpenBamFile(bam_file)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)
  chr_name <- names(bam_header[['targets']])
  chr_length <- unname(bam_header[['targets']])
  bam_header <- V2 <- V3 <- NULL
  test_list <- utils::read.csv(roi, sep = "\t", header = FALSE)

  mut_table <- function(V1){
    result <- which(chr_name == V1)
    return(result)
  }

  test <- unlist(sapply(test_list$V1, mut_table))
  test_list <- test_list %>% dplyr::mutate(V2 = V2 + 1, V3 = V3 + 1)
  test_list <- test_list %>%
    dplyr::mutate(chromosome = test) %>%
    dplyr::mutate(length = chr_length[chromosome])

  length <- test_list$length

  chrom_name <- test_list$V1
  reg_start <- test_list$V2
  reg_stop <- test_list$V3
  length <- test_list$length
  chromosome <- unlist(unname(test_list$chromosome))
  var_list <- list(chrom_name = chrom_name,
                   reg_start = reg_start,
                   reg_stop = reg_stop,
                   length = length,
                   chromosome = chromosome,
                   plot_output = plot_output,
                   path_to_RNAfold = path_to_RNAfold,
                   min_read_count = min_read_count,
                   genome= genome,
                   bam_file = bam_file,
                   roi = roi,
                   pi_pal = pi_pal,
                   si_pal = si_pal,
                   annotate_region = annotate_region,
                   weight_reads = weight_reads,
                   gtf_file = gtf_file,
                   write_fastas = write_fastas)
  return(var_list)
}
