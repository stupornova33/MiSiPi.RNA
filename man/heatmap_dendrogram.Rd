% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap_dendrogram.R
\name{heatmap_dendrogram}
\alias{heatmap_dendrogram}
\title{extracts read data and plots counts by read size
for multiple loci in a bed file}
\usage{
heatmap_dendrogram(bam_file, logfile, chrom_name, reg_start, reg_stop)
}
\arguments{
\item{bam_file}{a a string}

\item{logfile}{a string}

\item{chrom_name}{a list of chromosome names}

\item{reg_start}{a list of integers}

\item{reg_stop}{a list of integers}
}
\value{
a heatmap plot
}
\description{
extracts read data and plots counts by read size
for multiple loci in a bed file
}
