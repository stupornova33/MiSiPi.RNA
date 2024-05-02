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
#' @param weight_reads Determines whether read counts will be weighted. Valid options are "Top", "locus_norm", or "None". See MiSiPi documentation for descriptions of the weighting methods.
#' @param gtf_file a string corresponding to the path of genome annotation in 9-column GTF format.
#' @param write_fastas A string, "T" or "F". Optional. If "T", read pairs from functions will be written to file.
#' @param out_type The type of file for plots. Options are "png" or "pdf". Default is PDF.
#' @return a list
#' @export

set_vars <- function(roi, bam_file, genome, min_read_count, plot_output, path_to_RNAfold, pi_pal, si_pal, annotate_region, weight_reads, gtf_file = "F", write_fastas = "F", out_type = "pdf"){

  bam_obj <- OpenBamFile(bam_file)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)
  chr_name <- names(bam_header[['targets']])
  chr_length <- unname(bam_header[['targets']])
  bam_header <- V2 <- V3 <- NULL
  test_list <- utils::read.csv(roi, sep = "\t", header = FALSE)

  # get the working directory to write the logfile
  dir <- unlist(strsplit(roi, "\\/\\s*(?=[^\\/]+$)", perl=TRUE))[1]
  # assign indexes to the chromosomes names from the bed file
  # also checks to make sure the chromosome name from the bed file is actually in the genome
  # prints to a file

  res_list <- vector()
  na_idx <- vector()
  for(i in 1:nrow(test_list)){
    res <- which(chr_name == test_list$V1[i])
    if(identical(res, integer(0))){
      na_idx <- append(na_idx, i)
    } else {
      res_list <- append(res_list, res)
    }
  }

  if(length(na_idx) > 0){
    # remove any lines of bed file where chromosome was not in genome and print error to file.
    test_list <- test_list[-c(na_idx),]
    suppressWarnings(
      if(!file.exists("Error.log")){
        write(paste0("Chromosome at lines ", na_idx, " were not found in the genome. Please check.\n"), file = paste0(dir, "Error.log"), append = FALSE)
      } else {
        write(paste0("Chromosome at lines ", na_idx, " were not found in the genome. Please check.\n"), file = paste0(dir, "Error.log"), append = TRUE)
      }
    )
  }

  test_list <- test_list %>% dplyr::mutate(V2 = V2 + 1, V3 = V3 + 1)

  test_list <- test_list %>%
    dplyr::mutate(chromosome = res_list) %>%
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
                   write_fastas = write_fastas,
                   out_type = out_type)

  return(var_list)
}
