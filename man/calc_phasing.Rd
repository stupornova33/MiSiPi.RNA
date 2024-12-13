% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_phasing.R
\name{calc_phasing}
\alias{calc_phasing}
\title{calc phasing
plots output
takes two read dataframes
returns table of phasing scores}
\usage{
calc_phasing(df1, df2, n)
}
\arguments{
\item{df1}{a dataframe of summarized reads with a duplicates column}

\item{df2}{a dataframe of summarized reads with a duplicates column}

\item{n}{the number of bases to transform the read by}
}
\value{
phased_counts
}
\description{
calc phasing
plots output
takes two read dataframes
returns table of phasing scores
}
