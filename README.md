# MiSiPi-RNA
Characterization of small RNA pathways

#after installing devtools
#run devtools::install_github("stupornova33-MiSiPi-RNA")


library(misipi)

`%>%` <- magrittr::`%>%`
setwd("desired/working/directory")
#modify these variables


vars <- set_vars(bed_file = "path/to/bed", input_file = "path/to/bamfile", 
                 genome = "path/to/genome", min_read_count = 1, plot_output = "T", 
                 path_to_RNAfold = "path/to/ViennaRNA/RNAfold.exe", pi_pal = "BlYel", si_pal = "RdYlBl", annotate_bed = "T",
                 gff_file = "extdata/processed_dmel_annot.txt") #gff is included in misipi extdata



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
