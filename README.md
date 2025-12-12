<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/stupornova33/MiSiPi.RNA/assets/49455915/0f2caeae-12d2-4e87-9b85-44f568a1c44e">
  <source media="(prefers-color-scheme: light)" srcset="https://user-images.githubusercontent.com/63005660/236967995-82baabed-6ebf-45e1-a2d2-7e5ab27451a2.png">
  <img alt="MiSiPi R Package Logo" src="https://user-images.githubusercontent.com/63005660/236967995-82baabed-6ebf-45e1-a2d2-7e5ab27451a2.png">
</picture>

For more details about the package or to cite, please visit https://www.biorxiv.org/content/10.1101/2023.05.07.539760v1.

## MiSiPi.RNA
Characterization of small RNA pathways

### Installation and Basic Usage
You can find the full documentation and examples [here](https://github.com/stupornova33/MiSiPi.RNA/blob/main/documentation/Documentation.md).

In order to install MiSiPi.RNA, you must first install devtools and BiocManager:

```
install.packages("devtools")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  
devtools::install_github("stupornova33/MiSiPi.RNA")

library(MiSiPi.RNA)

```

### RNAfold
In order for this package to work, you must also have RNAfold and RNAplot (version >= 2.7.0) from the ViennaRNA
package installed. See https://www.tbi.univie.ac.at/RNA/ for installation.

#### Optional dependencies 
For converting the .ps/.eps output files from the siRNA/miRNA module to .png (see steps after "run_all" for examples), install [ImageMagick](https://imagemagick.org/index.php) and [ghostscript](https://www.ghostscript.com/).


### Input
The input for MiSiPi.RNA's main function is an object created by the set_vars() function. Running set_vars will always be the first step in using this package. Below is a description of each of the parameters that will be passed to set_vars(). These should be changed based on your needs.

- **roi**             - A bed file listing your regions of interest
- **bam_file**        - A BAM file of aligned reads. Index file must also be present
- **genome**          - A genome fasta file. Chromosome names must match the bed file
- **plot_output**     - (TRUE or FALSE) If TRUE, MiSiPi.RNA will output plots as pdfs
- **path_to_RNAfold** - Full path to RNAfold executable
- **path_to_RNAplot** - Full path to RNAplot executable
- **pi_pal**          - Palette option for the generated piRNA heatmap (see below)
- **si_pal**          - Palette option for the generated siRNA heatmap (see below)
- **annotate_region** - (TRUE or FALSE) Plots annotated gene features below the hairpin arc plot which is useful for characterizing cisNAT loci
- **weight_reads**    - Determines if read counts will be weighted. ("none", "locus_norm", "weight_by_prop", or an integer) 
- **gtf_file**        - Full path to a 9 column GTF file. Required only if annotate_region is TRUE. Note that the annotation column must contain features designated with "gene_id" or "transcript_id" (see https://www.ensembl.org/info/website/upload/gff.html?redirect=no#fields).
- **write_fastas**    - (TRUE or FALSE) If TRUE, MiSiPi.RNA will write read pairs from functions to a file. Default is FALSE
- **out_type**        - ("pdf" or "png") Specifies the output type. Default is "pdf"
- **use_bed_names**   - (Optional) Specifies whether column 4 of the region of interest BED file should be used as names for output files (Default = FALSE)
- **density_timeout** - (Optional) The density plot and calc_phasing functions can take awhile to run for highly abundant loci, and thus we set a default timeout of 3600 seconds (1 hour). If you wish to override this, provide a time in seconds. 

### Palettes:
Palette options are:
- "RdYlBl"
- "BlYel"
- "yelOrRed"
- "MagYel"
- "Greens"

Example usage:

```

library(MiSiPi.RNA)

vars <- set_vars(
    roi = "path/to/bed",
    bam_file = "path/to/bam", 
    genome = "path/to/genome",
    plot_output = TRUE, 
    path_to_RNAfold = "path/to/ViennaRNA/RNAfold.exe",
    path_to_RNAplot = "path/to/ViennaRNA/RNAplot.exe",
    pi_pal = "BlYel",
    si_pal = "RdYlBl",
    annotate_region = TRUE,
    weight_reads = "none",
    gtf_file = "path/to/gtf",
    write_fastas = FALSE,
    out_type = "pdf", 
    use_bed_names = FALSE,
    density_timeout = 3600
)



# run only for siRNAs
misipi_rna(vars, method = "siRNA")

```


The ```method``` parameter determines if your files will be processed for MicroRNA (```"miRNA"```), Piwi-interacting RNA (```"piRNA"```), Small interferring RNA (```"siRNA"```). By default, misipi_rna runs all three methods and combines the output plots. 

## Running with "all" method:
```
# run all methods
misipi_rna(vars)
```
In addition to processing files for miRNA, piRNA, and siRNA, the ```"all"``` method outputs a table with metrics and statistics which can be used for summarization or machine learning. See the [documentation](https://github.com/stupornova33/MiSiPi.RNA/blob/main/documentation/Documentation.html) for more details regarding values in table.

## Example Data Provided with MiSiPi.RNA
See examples/examples.Rmd for more.

```
library(MiSiPi.RNA)

roi <- system.file("extdata", "dmel_mir-bantam.bed", package = "MiSiPi.RNA")
bam_file <- system.file("extdata", "dmel_bantam.bam", package = "MiSiPi.RNA")
genome <- system.file("extdata", "dmel_chr3L.fasta", package = "MiSiPi.RNA")
annot <- system.file("extdata", "dmel_example.gtf", package = "MiSiPi.RNA")
path_to_RNAfold = "" # Replace with your path to RNAfold and RNAplot here
path_to_RNAplot = ""

vars <- set_vars(
  roi = roi,
  bam_file = bam_file,
  genome = genome,
  path_to_RNAfold = path_to_RNAfold,
  path_to_RNAplot = path_to_RNAplot,  
  pi_pal = "BlYel",
  si_pal = "RdYlBl",
  annotate_region = TRUE,
  gtf_file = annot,
  out_type = "png")

misipi_rna(vars, method = "miRNA")

```
## To convert .ps/.eps files to .png 

```
eps2png(path_to_magick_exe, file_dir)
```
where path_to_magick_exe is the full path to the binary executable and file_dir is the folder containing the files to convert. This will also be the output folder.

## Use built-in machine learning model to characterize loci:
See full [documentation](https://github.com/stupornova33/MiSiPi.RNA/blob/main/documentation/Documentation.html) for more details. 

```
ml_probability(path_to_table = "full/path/to/run_all/", table = "table_ml.txt")
```
The ```path_to_table``` parameter is the path to the directory that was created when ```misipi_rna()``` was run with the default method. It will be called "run_all/".


## To make a plot summary and sortable table of results:

```
make_html_summary(path_to_tables = "full/path/to/run_all/", type = c("siRNA", "piRNA" or "miRNA"), ml_plots = FALSE)
```
The ```path_to_tables``` parameter is the path to the directory that was created when ```misipi_rna()``` was run with the default method. It will be called "run_all/".

The ```ml_plots``` parameter is intended for users that have already run ```ml_probability()```. The default for this option is false. 





