#' plot miRNA coverage
#' takes a data frame of coverage over a sequence of nucleotides

#' @param bed_file a tab-delineated file with three columns: chrom, start, stop
#' @param input_file a BAM file
#' @param genome_file a fasta file with chromosome sequences
#' @param plot_output an optional string. Default = TRUE
#' @param path_to_RNAfold a string
#' @importFrom Rcpp sourceCpp
#' @return plots

#' @export

run_miRNA_function <- function(bed_file, input_file, genome_file, plot_output = plot_output, path_to_RNAfold){
   `%>%` <- magrittr::`%>%`
   
   logfile = "miRNA_logfile.txt"

   dir <- 'miRNA_outputs/'
   if(!dir.exists(dir) == TRUE) dir.create(dir)
   
   if(!file.exists(logfile) == TRUE) file.create(paste0(dir, logfile))

   mapply(miRNA_function, vars[[1]], vars[[2]], vars[[3]], vars[[5]], vars[[4]], "+", 1, vars[[9]], vars[[10]], logfile, dir, vars[[6]], vars[[7]])
   mapply(miRNA_function, vars[[1]], vars[[2]], vars[[3]], vars[[5]], vars[[4]], "-", 1, vars[[9]], vars[[10]], logfile, dir, vars[[6]], vars[[7]])
   
   
}