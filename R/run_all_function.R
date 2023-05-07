#' run the run_all function
#' @param vars a list
#' @importFrom Rcpp sourceCpp
#' @return plots

#' @export

run_all_function <- function(vars){
   chrom_name <- vars[[1]]
   reg_start <- vars[[2]]
   reg_stop <- vars[[3]]
   chromosome <- vars[[5]]
   length <- vars[[4]]
   input_file <- vars[[10]]
   genome_file <- vars[[9]]
   min_read_count <- vars[[8]]
   si_pal <- vars[[13]]
   pi_pal <- vars[[12]]
   plot_output <- vars[[6]]
   path_to_RNAfold <- vars[[7]]
   bed_file <- vars[[11]]
   
   `%>%` <- magrittr::`%>%`
   #bam_obj <- OpenBamFile(input_file)
   #bam_header <- Rsamtools::scanBamHeader(bam_obj)
   #chr_name <- names(bam_header[['targets']])
   #chr_length <- unname(bam_header[['targets']])
   
   #bam_header <- NULL
   #test_list <- read.csv(bed_file, sep = "\t", header = FALSE)
   
   #mut_table <- function(V1){
   #   result <- which(chr_name == V1)
      
   #}
   
   #test <- unlist(sapply(test_list$V1, mut_table))
   #test_list <- test_list %>% dplyr::mutate(V2 = V2 + 1, V3 = V3 + 1)
   #test_list <- test_list %>%
   #   dplyr::mutate(chromosome = test) %>%
   #   dplyr::mutate(length = chr_length[chromosome])
   

   mapply(new_run_all, chrom_name, reg_start, reg_stop, chromosome, length, input_file, genome_file, min_read_count, 
          si_pal, pi_pal, plot_output, path_to_RNAfold, bed_file)
   
}

