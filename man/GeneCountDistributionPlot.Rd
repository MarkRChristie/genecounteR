% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneCountDistributionPlot.R
\name{GeneCountDistributionPlot}
\alias{GeneCountDistributionPlot}
\title{Gene Count Distribution Plots}
\usage{
GeneCountDistributionPlot(dat, threshold)
}
\arguments{
\item{dat}{gene count matrix}

\item{threshold}{percent of data to retain (e.g., if threshold set to 0.15, top 15\% (mean expression) of data retained); 0 < threshold < 1}
}
\value{
A filtered count matrix
}
\description{
Retain top x percent of data; results illustrated on plot (red dashed line illustrates threshold)
}
