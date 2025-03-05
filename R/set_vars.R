#' Function to process and set the arguments passed by the user for each function
#' @param roi The path to a BED file of loci of interest
#' @param bam_file The path to a BAM file
#' @param genome The path to a genome Fasta file
#' @param path_to_RNAfold The full path to the RNAfold binary executable.
#' @param path_to_RNAplot The full path to the RNAplot binary executable.
#' @param plot_output Determines whether the program will output plots as PDFs. Expected input is TRUE or FALSE.
#' @param pi_pal The color palette to use for the piRNA heatmap plot. Valid options are "RdYlBl", "BlYel", "yelOrRed", "MagYel", and "Greens".
#' @param si_pal The color palette to use for the siRNA heatmap plot. Valid options are "RdYlBl", "BlYel", "yelOrRed", "MagYel", and "Greens".
#' @param weight_reads Determines whether read counts will be weighted. Valid options are "Top", "locus_norm", or "None". See MiSiPi documentation for descriptions of the weighting methods.
#' @param write_fastas TRUE or FALSE. Optional. If TRUE, read pairs from functions will be written to file.
#' @param annotate_region Determines whether the program will plot genomic features of interest found in the GTF annotation file. If TRUE, a GTF file must be provided as the "gtf_file" argument.
#' @param gtf_file a string corresponding to the path of genome annotation in 9-column GTF format. Default is FALSE unless annotate_regions == TRUE.
#' @param out_type The type of file for plots. Options are "png" or "pdf". Default is PDF.
#' @return a list
#' @export

set_vars <- function(roi, bam_file, genome, 
                     path_to_RNAfold, path_to_RNAplot, plot_output = TRUE, 
                     pi_pal = c("RdYlBl", "BlYel", "yelOrRed", "MagYel", "Greens"),
                     si_pal = c("RdYlBl", "BlYel", "yelOrRed", "MagYel", "Greens"),
                     weight_reads = c("None", "top", "locus_norm", "none", "Top", "Locus_Norm"), 
                     write_fastas = FALSE, annotate_region = FALSE, gtf_file = FALSE,
                     out_type = c("pdf", "png", "PDF", "PNG")) {
  ## Parameter Validation
  # roi
  stopifnot("Parameter `roi` must be a valid filepath to a BED file." = file.exists(roi))
  bed_columns_vector <- utils::count.fields(roi, sep = "\t")
  stopifnot("Bed file (roi) must have the same number of columns in each line." = length(unique(bed_columns_vector)) == 1)
  number_of_bed_columns <- bed_columns_vector[1]
  stopifnot("Bed file (roi) must have 3 columns and be tab separated." = number_of_bed_columns >= 3)
  # bam_file
  stopifnot("Parameter `bam_file` must have a .bam extension." = tools::file_ext(bam_file) == "bam")
  stopifnot("Parameter `bam_file` must be a valid filepath to a BAM file." = file.exists(bam_file))
  bai_file <- paste0(bam_file, ".bai")
  stopifnot("A corresponding .bai index file must be present in the same directory as the .bam file" = file.exists(bai_file))
  # genome
  stopifnot("Parameter `genome` must be a valid filepath to a genome Fasta file." = file.exists(genome))
  # plot_output
  stopifnot("Parameter `plot_output` only accepts TRUE or FALSE." = is.logical(plot_output))
  # path_to_RNAfold
  stopifnot("Parameter `path_to_RNAfold` must be a valid filepath to RNAfold." = file.exists(path_to_RNAfold))
  stopifnot("Parameter `path_to_RNAplot` must be a valid filepath to RNAplot." = file.exists(path_to_RNAplot))
  # pi_pal - Dependent on plot_output
  pi_pal <- match.arg(pi_pal)
  # si_pal - Dependent on plot_output
  si_pal <- match.arg(si_pal)
  # annotate_region
  stopifnot("Parameter `annotate_region` only accepts TRUE or FALSE." = is.logical(annotate_region))
  # weight_reads
  weight_reads <- match.arg(weight_reads)
  # gtf_file - Dependent on annotate_region
  if (annotate_region == TRUE) {
    stopifnot("Parameter `gtf_file` must be provided when `annotate_region` is TRUE." = !missing(gtf_file))
    stopifnot("Parameter `gtf_file` must be a valid filepath to a 9 column gtf file." = file.exists(gtf_file))
    gtf_columns_vector <- utils::count.fields(gtf_file, sep = "\t")
    stopifnot("gtf_file must have the same number of columns in each line." = length(unique(gtf_columns_vector)) == 1)
    # waiting to include the number of column check
    # number_of_gtf_columns <- gtf_columns_vector[1]
    # stopifnot("gtf_file must have 9 columns and be tab separated." = number_of_gtf_columns == 9)
  }
  # write_fastas
  stopifnot("Parameter `write_fastas` only accepts TRUE or FALSE." = is.logical(write_fastas))
  # out_type
  out_type <- match.arg(out_type)
  out_type <- tolower(out_type)
  ## End Parameter Validation

  bam_obj <- .open_bam(bam_file)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)
  chr_name <- names(bam_header[["targets"]])
  chr_length <- unname(bam_header[["targets"]])
  bam_header <- V2 <- V3 <- NULL
  test_list <- utils::read.csv(roi, sep = "\t", header = FALSE)

  # get the working directory to write the logfile
  dir <- unlist(strsplit(roi, "\\/\\s*(?=[^\\/]+$)", perl = TRUE))[1]
  # assign indexes to the chromosomes names from the bed file
  # also checks to make sure the chromosome name from the bed file is actually in the genome
  # prints to a file

  res_list <- vector()
  na_idx <- vector()
  for (i in 1:nrow(test_list)) {
    res <- which(chr_name == test_list$V1[i])
    if (identical(res, integer(0))) {
      na_idx <- append(na_idx, i)
    } else {
      res_list <- append(res_list, res)
    }
  }

  if (length(na_idx) > 0) {
    # remove any lines of bed file where chromosome was not in genome and print error to file.
    test_list <- test_list[-c(na_idx), ]
    suppressWarnings(
      if (!file.exists("Error.log")) {
        write(paste0("Chromosome at lines ", na_idx, " were not found in the genome. Please check.\n"), file = paste0(dir, "Error.log"), append = FALSE)
      } else {
        write(paste0("Chromosome at lines ", na_idx, " were not found in the genome. Please check.\n"), file = paste0(dir, "Error.log"), append = TRUE)
      }
    )
  }

  # Convert the bed file coordinates to 1 based for compatibility with other tools
  # This will be converted back using revert_positions() when writing results that reference the coordinates
  test_list <- test_list %>%
    dplyr::mutate(
      V2 = V2 + 1,
      V3 = V3 + 1
    )

  test_list <- test_list %>%
    dplyr::mutate(chromosome = res_list) %>%
    dplyr::mutate(length = chr_length[chromosome])

  length <- test_list$length

  chrom_name <- test_list$V1
  reg_start <- test_list$V2
  reg_stop <- test_list$V3
  length <- test_list$length
  chromosome <- unlist(unname(test_list$chromosome))
  var_list <- list(
    chrom_name = chrom_name,
    reg_start = reg_start,
    reg_stop = reg_stop,
    length = length,
    chromosome = chromosome,
    plot_output = plot_output,
    path_to_RNAfold = path_to_RNAfold,
    path_to_RNAplot = path_to_RNAplot,
    genome = genome,
    bam_file = bam_file,
    roi = roi,
    pi_pal = pi_pal,
    si_pal = si_pal,
    annotate_region = annotate_region,
    weight_reads = weight_reads,
    gtf_file = gtf_file,
    write_fastas = write_fastas,
    out_type = out_type
  )

  return(var_list)
}
