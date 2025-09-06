#' Function to process and set the arguments passed by the user for each function
#' @param roi The path to a BED file of loci of interest
#' @param bam_file The path to a BAM file
#' @param genome The path to a genome Fasta file
#' @param path_to_RNAfold The full path to the RNAfold binary executable.
#' @param path_to_RNAplot The full path to the RNAplot binary executable.
#' @param plot_output Determines whether the program will output plots as PDFs.
#'   Expected input is TRUE or FALSE.
#' @param pi_pal The color palette to use for the piRNA heatmap plot.
#'   Valid options are "RdYlBl", "BlYel", "yelOrRed", "MagYel", and "Greens".
#' @param si_pal The color palette to use for the siRNA heatmap plot.
#'   Valid options are "RdYlBl", "BlYel", "yelOrRed", "MagYel", and "Greens".
#' @param weight_reads Determines whether read counts will be weighted.
#'   Valid options are "Top", "locus_norm", or "None".
#'   See MiSiPi documentation for descriptions of the weighting methods.
#' @param write_fastas TRUE or FALSE. Optional.
#'   If TRUE, read pairs from functions will be written to file.
#' @param annotate_region Determines whether the program will plot genomic
#'   features of interest found in the GTF annotation file.
#'   If TRUE, a GTF file must be provided as the "gtf_file" argument.
#' @param gtf_file a string corresponding to the path of genome
#'   annotation in 9-column GTF format.
#'   Default is FALSE unless annotate_regions == TRUE.
#' @param out_type The type of file for plots.
#'   Options are "png" or "pdf". Default is PDF.
#' @param use_bed_names a boolean indicating if the name column in
#'   the bed file should be used for results.
#'   If TRUE, bed line names will be used unless column 4 is not present in the bed file.
#'   In that case, this parameter will toggle to FALSE.
#'   If FALSE, results will be named using a region string in the format: chr-start_stop
#' @return a list
#' @export

set_vars <- function(roi, bam_file, genome, 
                     path_to_RNAfold, path_to_RNAplot, plot_output = TRUE, 
                     pi_pal = c("RdYlBl", "BlYel", "yelOrRed", "MagYel", "Greens"),
                     si_pal = c("RdYlBl", "BlYel", "yelOrRed", "MagYel", "Greens"),
                     weight_reads = c("None", "top", "locus_norm", "none", "Top", "Locus_Norm"), 
                     write_fastas = FALSE, annotate_region = FALSE, gtf_file = FALSE,
                     out_type = c("pdf", "png", "PDF", "PNG"),
                     use_bed_names = FALSE) {
  #### Parameter Validation ####
  # roi - bed file
  stopifnot("Parameter `roi` must be a valid filepath to a BED file." = file.exists(roi))
  bed_columns_vector <- utils::count.fields(roi, sep = "\t")
  stopifnot("Bed file (roi) must have the same number of columns in each line." = length(unique(bed_columns_vector)) == 1)
  number_of_bed_columns <- bed_columns_vector[1]
  stopifnot("Bed file (roi) must have 3 columns and be tab separated." = number_of_bed_columns >= 3)
  
  # If there are more than 3 columns, assume the 4th column is region name
  bed_names_present <- ifelse(number_of_bed_columns > 3, TRUE, FALSE)
  
  # If use_bed_names has been set to TRUE, ensure that a 4th column is actually present in the bed file
  if (use_bed_names & !bed_names_present) {
    cli::cli_warn("Bed file only has 3 columns, so bed file names will be replaced with [chr-start_stop] in results.")
    use_bed_names <- FALSE
  }
  
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
    gtf_columns_vector <- utils::count.fields(gtf_file, sep = "\t", quote = "")
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
  
  #### Vienna Software Version Validation ####
  os <- Sys.info()["sysname"]
  
  # To be used for RNAFold and RNAPlot version checking
  getVersion <- function(os, supplied_path) {
    if (os == "Windows") {
      rnafold_version <- system2(
        command = supplied_path,
        args = "--version",
        stdout = TRUE,
        wait = TRUE,
        invisible = TRUE
      )
    } else if (os == "Linux" | os == "Darwin") {
      rnafold_version <- system(
        command = paste(supplied_path, "--version", sep = " "),
        intern = TRUE
      )
    } else {
      warning("Operating system is not Windows or Linux/Darwin. Report as a bug if you think this is a mistake.")
      return(NULL)
    }
    return(rnafold_version)
  }
  
  # Attempt a system call to get RNAFold version
  tryCatch(
    rnafold_version <- getVersion(os, path_to_RNAfold),
    error = function(e) {
      if (grepl("not found", e$message)) {
        msg <- paste(cli::col_br_red("PATH ERROR:"), path_to_RNAfold, "not found. Provide correct path to RNAfold.")
      } else {
        msg <- e$message
      }
      cli::cli_abort(msg, call = NULL)
    }
  )
  
  # Attempt a system call to get RNAPlot version
  tryCatch(
    rnaplot_version <- getVersion(os, path_to_RNAplot),
    error = function(e) {
      if (grepl("not found", e$message)) {
        msg <- paste(cli::col_br_red("PATH ERROR:"), path_to_RNAplot, "not found. Provide correct path to RNAplot.")
      } else {
        msg <- e$message
      }
      cli::cli_abort(msg, call = NULL)
    }
  )
  
  # Validate RNAFold
  rnafold_version <- stringr::str_split_1(rnafold_version, " ")[2]
  # str_split uses regex to split, so when splitting by ".", brackets must be placed around it
  version_split <- stringr::str_split_1(rnafold_version, "[.]")
  major <- as.integer(version_split[1])
  minor <- as.integer(version_split[2])
  
  if (major < 2L || (major >= 2L && minor < 7L)) {
    msg <- paste(cli::col_br_red("VERSION ERROR"), "RNAFold version:", rnafold_version, "---- Must be at least 2.7.0")
    cli::cli_abort(msg, call = NULL)
  }
  
  # Validate RNAPlot
  rnaplot_version <- stringr::str_split_1(rnaplot_version, " ")[2]
  # str_split uses regex to split, so when splitting by ".", brackets must be placed around it
  version_split <- stringr::str_split_1(rnaplot_version, "[.]")
  major <- as.integer(version_split[1])
  minor <- as.integer(version_split[2])
  
  if (major < 2L || (major >= 2L && minor < 7L)) {
    msg <- paste(cli::col_br_red("VERSION ERROR"), "RNAPlot version:", rnaplot_version, "---- Must be at least 2.7.0")
    cli::cli_abort(msg, call = NULL)
  }
  
  # Scan bam file header for chromosome names
  bam_obj <- .open_bam(bam_file)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)
  chr_name <- names(bam_header[["targets"]])
  bam_header <- V2 <- V3 <- NULL
  .close_bam(bam_obj)
  
  # Read in bed file
  bed_lines <- utils::read.csv(roi, sep = "\t", header = FALSE)

  # assign indexes to the chromosomes names from the bed file
  # also checks to make sure the chromosome name from the bed file is actually in the bam file
  res_list <- vector()
  na_idx <- vector()
  
  for (i in 1:nrow(bed_lines)) {
    res <- which(chr_name == bed_lines$V1[i])
    if (identical(res, integer(0))) {
      na_idx <- append(na_idx, i)
    } else {
      res_list <- append(res_list, res)
    }
  }

  stopifnot("There are no matching chromosomes between bed file and bam file." = length(res_list) > 0)
  
  if (length(na_idx) > 0) {
    # remove any lines of bed file where chromosome was not in genome and print error to file.
    bed_lines <- bed_lines[-c(na_idx), ]
    
    warning_message <- paste("Chromosome at lines", na_idx, "were not found in the genome.\n")
    cli::cli_warn(warning_message)
  }
  
  res_list <- na_idx <- NULL

  # Convert the bed file coordinates to 1 based for compatibility with Rsamtools
  # Bed files use zero-based half open coordinates, so only the start position needs to be incremented
  # Coordinates will be reverted back to the original bed file coordinates when writing output
  bed_lines <- bed_lines %>%
    dplyr::mutate(
      # Increment start position by 1
      V2 = V2 + 1,
      V3 = V3
    )

  chrom_name <- bed_lines$V1
  reg_start <- bed_lines$V2
  reg_stop <- bed_lines$V3
  
  if (use_bed_names) {
    iteration_output <- prefix <- bed_lines$V4
  } else {
    iteration_output <- .get_bed_region_string(chrom_name, reg_start, reg_stop)
    prefix <- .get_region_string(chrom_name, reg_start, reg_stop)
  }
  
  bed_lines <- NULL
  
  var_list <- list(
    chrom_name = chrom_name,
    reg_start = reg_start,
    reg_stop = reg_stop,
    prefix = prefix,
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
    out_type = out_type,
    use_bed_names = use_bed_names,
    iteration_output = iteration_output
  )

  return(var_list)
}
