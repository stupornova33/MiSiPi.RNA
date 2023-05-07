

#Rcpp::Rcpp.package.skeleton("misipi",cpp_files = "C:/Users/tmjar/Desktop/misipi/src/cpp_functions.cpp", force = TRUE)
#restart R session
setwd("path/to/misipi/folder")
devtools::install("misipi", force = TRUE)

library(misipi)

`%>%` <- magrittr::`%>%`
setwd("desired/working/directory")
vars <- set_vars(bed_file = "path/to/bed", input_file = "path/to/bamfile", 
                 genome = "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna", min_read_count = 1, plot_output = "T", 
                 path_to_RNAfold = "path/to/ViennaRNA/RNAfold.exe", pi_pal = "BlYel", si_pal = "RdYlBl", annotate_bed = "T",
                 gff_file = "processed_dmel_annot.txt")


#to run one line of bed file for testing
chrom_name <- vars[[1]][1]
reg_start <- vars[[2]][1]
reg_stop <- vars[[3]][1]
length <- vars[[4]][1]
input_file <- vars[[10]]
genome_file <- vars[[9]]
min_read_count <- 1
path_to_RNAfold <- vars[[7]]
plot_output <- "T"
pal <- "RdYlBl"
dir <- "C:/Users/path/to/dir"
logfile <- 'logfile.txt'
bed_file <- vars[[11]]
annotate_bed <- "F"
gff_file <- vars[[15]]

#to run all lines of bed file
run_miRNA_function(vars)

run_piRNA_function(vars)


run_siRNA_function(vars)


run_phased_piRNA(vars)


