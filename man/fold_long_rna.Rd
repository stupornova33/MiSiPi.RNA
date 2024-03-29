% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fold_long_rna.R
\name{fold_long_rna}
\alias{fold_long_rna}
\title{function to run RNAfold
processes output of RNA fold to get MFE and vien struct
returns list of values for each region}
\usage{
fold_long_rna(chrom_name, start, stop, converted, path_to_RNAfold)
}
\arguments{
\item{chrom_name}{a string}

\item{start}{a whole number}

\item{stop}{a whole number}

\item{converted}{a vector containing a sequence}

\item{path_to_RNAfold}{a string}
}
\value{
list
}
\description{
function to run RNAfold
processes output of RNA fold to get MFE and vien struct
returns list of values for each region
}
