% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{get_phased_dist}
\alias{get_phased_dist}
\title{get_phased_dist}
\usage{
get_phased_dist(startX, endX, startY, endY, dfsize)
}
\arguments{
\item{startX}{A vector of ints}

\item{endX}{A vector of ints}

\item{startY}{A vector of ints}

\item{endY}{A vector of ints}

\item{dfsize}{An integer}
}
\value{
df A data.frame containing the start positions of both reads, their widths, and the distance
}
\description{
This function takes in vectors of start and stop positions for read groups 1 and 2
It calculates the distance between the end of one group and start of another and returns
The results in a data.frame
}
