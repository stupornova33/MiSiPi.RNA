% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_sizes_only.R
\name{plot_sizes_only}
\alias{plot_sizes_only}
\title{plot_sizes_only
A function which will process the reads at a locus but only plot the size distribution and will not run the miRNA, siRNA or piRNA modules.}
\usage{
plot_sizes_only(wkdir, geno, bam_file, bed_file)
}
\arguments{
\item{wkdir}{A string - the path of the output directory.}

\item{geno}{A string specifying the genotype of the sample, for appending as a suffix to the output file.}

\item{bam_file}{A string - the name of the BAM file.}

\item{bed_file}{A string - the name of the 3-column BED file.}
}
\value{
results
}
\description{
plot_sizes_only
A function which will process the reads at a locus but only plot the size distribution and will not run the miRNA, siRNA or piRNA modules.
}
