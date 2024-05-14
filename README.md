<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/stupornova33/MiSiPi.RNA/assets/49455915/0f2caeae-12d2-4e87-9b85-44f568a1c44e">
  <source media="(prefers-color-scheme: light)" srcset="https://user-images.githubusercontent.com/63005660/236967995-82baabed-6ebf-45e1-a2d2-7e5ab27451a2.png">
  <img alt="MiSiPi R Package Logo" src="https://user-images.githubusercontent.com/63005660/236967995-82baabed-6ebf-45e1-a2d2-7e5ab27451a2.png">
</picture>

For more details about the package or to cite, please visit https://www.biorxiv.org/content/10.1101/2023.05.07.539760v1.

## MiSiPi.RNA
Characterization of small RNA pathways

### Installation
In order to install MiSiPi.RNA, you must first install devtools and BiocManager:

```
install.packages("devtools")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  
devtools::install_github("stupornova33/MiSiPi.RNA")

library(MiSiPi.RNA)

```

### RNAfold
In order for this package to work, you must also have RNAfold from the ViennaRNA
package installed. You will need the path to the RNAfold executable.
See https://www.tbi.univie.ac.at/RNA/ for installation.

### Input
The input for any of MiSiPi.RNA's main functions is an object created by the
set_vars() function. Running set_vars will always be the first step in using
this package. Below is a description of each of the parameters that will be
passed to set_vars(). These should be changed based on your needs.

- **roi**             - A bed file listing your regions of interest
- **bam_file**        - A BAM file of aligned reads. Index file must also be present
- **genome**          - A genome fasta file. Chromosome names must match the bed file
- **min_read_count**  - This filters out loci with low mapping reads. Defaults to 1
- **plot_output**     - (TRUE or FALSE) If TRUE, MiSiPi.RNA will output plots as pdfs
- **path_to_RNAfold** - Full path to RNAfold executable
- **pi_pal**          - Palette option for the generated piRNA heatmap (see below)
- **si_pal**          - Palette option for the generated siRNA heatmap (see below)
- **annotate_region** - (TRUE or FALSE) Plots annotated gene features below the hairpin arc plot which is useful for characterizing cisNAT loci
- **weight_reads**    - Determines if read counts will be weighted. ("Top", "locus_norm", or "None") 
- **gtf_file**        - Full path to a 9 column GTF file. Required only if annotate_region is TRUE
- **write_fastas**    - (TRUE or FALSE) If TRUE, MiSiPi.RNA will write read pairs from functions to a file. Default is FALSE
- **out_type**        - ("pdf" or "png") Specifies the output type. Default is "pdf"


```
vars <- set_vars(roi = "path/to/bed",
                bam_file = "path/to/bam", 
                genome = "path/to/genome",
                min_read_count = 1,
                plot_output = TRUE, 
                path_to_RNAfold = "path/to/ViennaRNA/RNAfold.exe",
                pi_pal = "BlYel",
                si_pal = "RdYlBl",
                annotate_region = TRUE,
                weight_reads = "None",
                gtf_file = "path/to/gtf",
                write_fastas = FALSE,
                out_type = "pdf")

```

### Palettes:
Palette options are:
- "RdYlBl"
- "BlYel"
- "yelOrRed"
- "MagYel"
- "Greens"


## To run all lines of bed file:

```
miRNA_function(vars)


piRNA_function(vars)


siRNA_function(vars)

# To run the above functions all at once
misipi_rna(vars)

```
