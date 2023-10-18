
![RpackageLogo](https://user-images.githubusercontent.com/63005660/236967995-82baabed-6ebf-45e1-a2d2-7e5ab27451a2.png)

For more details about the package or to cite, please visit https://www.biorxiv.org/content/10.1101/2023.05.07.539760v1.

## MiSiPi.RNA
Characterization of small RNA pathways

after installing devtools and BiocManager run 

```
devtools::install_github("stupornova33/MiSiPi.RNA")

library(MiSiPi.RNA)

`%>%` <- magrittr::`%>%`

setwd("desired/working/directory")
```

## Modify these variables.
Bed_file is your list of regions of interest in BED format. Input file is the BAM file of aligned reads (must be indexed). Path to genome is the path to your genome fasta file (this must have the same chromosome names that are found in the BED file). Min read count filters out loci with low reads mapping; default is 1. 
RNA fold must be installed on your system, but path_to_RNAfold is only the path to the RNAfold binary .exe.

```
ROI <- set_vars(roi = "path/to/bed", bam_file = "path/to/bamfile", 
                genome = "path/to/genome", min_read_count = 1, plot_output = "T", 
                path_to_RNAfold = "path/to/ViennaRNA/RNAfold.exe", pi_pal = "BlYel", si_pal = "RdYlBl", annotate_bed = "T",
                bed_file = "extdata/processed_dmel_annot.txt", weight_reads = "F") #sample BED is included in misipi extdata

```

### Palettes:
Pi_pal is the piRNA heatmap plot, Si_pal is the siRNA heatmap plot. 
Palette options are: "RdYlBl", "BlYel", "yelOrRed", "MagYel", and "Greens". 

### Locus annotation: 
If annotate_bed = "T", a BED file must be supplied to the bed_file argument. This will plot annotated gene features below the hairpin arc plot which is useful for characterizing cisNAT loci. 


## to run all lines of bed file

```
miRNA_function(vars)


piRNA_function(vars)


siRNA_function(vars)

```
