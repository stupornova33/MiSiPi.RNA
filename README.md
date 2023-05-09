
![RpackageLogo](https://user-images.githubusercontent.com/63005660/236967995-82baabed-6ebf-45e1-a2d2-7e5ab27451a2.png)



## MiSiPi-RNA
Characterization of small RNA pathways

after installing devtools run 
>devtools::install_github("stupornova33-MiSiPi-RNA")


>library(misipi)

>`%>%` <- magrittr::`%>%`
>setwd("desired/working/directory")

modify these variables


>vars <- set_vars(bed_file = "path/to/bed", input_file = "path/to/bamfile", 
>                 genome = "path/to/genome", min_read_count = 1, plot_output = "T", 
>                 path_to_RNAfold = "path/to/ViennaRNA/RNAfold.exe", pi_pal = "BlYel", si_pal = "RdYlBl", annotate_bed = "T",
>                 gff_file = "extdata/processed_dmel_annot.txt") #gff is included in misipi extdata


##to run all lines of bed file


>run_miRNA_function(vars)


>run_piRNA_function(vars)


>run_siRNA_function(vars)


>run_phased_piRNA(vars)
