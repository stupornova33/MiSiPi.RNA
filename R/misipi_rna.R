#' Wrapper function that calls core functions
#' @param vars a list of parameters created by set_vars
#' @param method a string indicating which processing method should be used
#' @return plots

#' @export

misipi_rna <- function(vars, method = c("all", "miRNA", "piRNA", "siRNA"), outdir_name) {
  method <- match.arg(method)
  
  # Create base output directory if it doesn't already exist
  # Directory will be made in directory from which misipi was called
  OUTPUT_BASE_DIR <- "misipi_output"
  if (!dir.exists(OUTPUT_BASE_DIR)) {
    dir.create(OUTPUT_BASE_DIR)
  }
  
  # Create current run's output directory inside the base output directory
  now <- format(lubridate::now(), "%m-%d-%Y_%H%M%S")
  
  output_dir <- file.path(OUTPUT_BASE_DIR, outdir_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Create a more human readable timestamp for the metadata file
  now_hr <- gsub("_", " ", now)
  now_hr <- sub("(\\d{2})(\\d{2})(\\d{2})$", "\\1:\\2:\\3", now_hr) 
  
  # Create a metadata object to track how the program was run
  metadata <- data.frame(Data = character(20))
  
  metadata$Data[1] <- paste0("Version: ", toString(packageVersion("MiSiPi.RNA")))
  metadata$Data[2] <- paste0("Started At: ", now_hr)
  metadata$Data[3] <- paste0("Output Directory: ", file.path(getwd(), output_dir))
  metadata$Data[4] <- paste0("Method: ", method)
  metadata$Data[5] <- paste0("Bam File: ", vars$bam_file)
  metadata$Data[6] <- paste0("Genome File: ", vars$genome)
  metadata$Data[7] <- paste0("Bed File: ", vars$roi)
  metadata$Data[8] <- paste0("Using Bed File Names: ", vars$use_bed_names)
  metadata$Data[9] <- paste0("GTF File: ", vars$gtf_file)
  metadata$Data[10] <- paste0("RNAfold Path: ", vars$path_to_RNAfold)
  metadata$Data[11] <- paste0("RNAplot Path: ", vars$path_to_RNAplot)
  metadata$Data[12] <- paste0("Weighting: ", vars$weight_reads)
  metadata$Data[13] <- paste0("Annotated: ", vars$annotate_region)
  metadata$Data[14] <- paste0("Plotted: ", vars$plot_output)
  metadata$Data[15] <- paste0("Plot Type: ", vars$out_type)
  metadata$Data[16] <- paste0("Fastas Written: ", vars$write_fastas)
  metadata$Data[17] <- paste0("piRNA Palette: ", vars$pi_pal)
  metadata$Data[18] <- paste0("siRNA Palette: ", vars$si_pal)
  metadata$Data[19] <- paste0("Read Density Timeout: ", vars$density_timeout)
  metadata$Data[20] <- paste0("Interactive: ", interactive())
  
  metadata_file <- file.path(output_dir, "run_info.txt")
  
  write.table(metadata, metadata_file, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  if (method == "all") {
    run_all(vars, output_dir)
  } else if (method == "miRNA") {
    miRNA(vars, output_dir)
  } else if (method == "piRNA") {
    piRNA(vars, output_dir)
  } else if (method == "siRNA") {
    siRNA(vars, output_dir)
  } else {
    stop(paste("`method`:", method, "is not valid."))
  }
}
