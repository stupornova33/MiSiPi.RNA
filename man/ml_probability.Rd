% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ml_probability.R
\name{ml_probability}
\alias{ml_probability}
\title{ml_probability}
\usage{
ml_probability(path_to_table, table, library_size)
}
\arguments{
\item{path_to_table}{A string specifying the full path to the folder which contains the _ml table. This directory will be the location of the probability plots.}

\item{table}{A string specifying the full path to the _ml table produced by misipi_rna}

\item{library_size}{An integer specifying the number of mapped reads from the BAM file that was used to create the table. This is used for RPM-mapped count normalization.}
}
\value{
nothing
}
\description{
ml_probability
}
