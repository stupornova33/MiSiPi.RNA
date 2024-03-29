% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_hp_overlaps.R
\name{find_hp_overlaps}
\alias{find_hp_overlaps}
\title{function to find hairpin overlaps
takes two data tables of reads with start and stop positions
converts to gr objects and finds overlaps, transforms positions back to original
returns overlap table}
\usage{
find_hp_overlaps(r1_dt, r2_dt, n)
}
\arguments{
\item{r1_dt}{a data frame containing chrom name, start, and stop}

\item{r2_dt}{a data frame containing chrom name, start, and stop}

\item{n}{the number of bases to transform the read by}
}
\value{
overlap table
}
\description{
function to find hairpin overlaps
takes two data tables of reads with start and stop positions
converts to gr objects and finds overlaps, transforms positions back to original
returns overlap table
}
