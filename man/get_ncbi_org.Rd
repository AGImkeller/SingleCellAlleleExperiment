% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SingleCellAlleleExperiment.R
\name{get_ncbi_org}
\alias{get_ncbi_org}
\title{Get NCBI genes using the org.HS.db package}
\usage{
get_ncbi_org(sce)
}
\arguments{
\item{sce}{A \code{\link{SingleCellExperiment}} object.}
}
\value{
A list of character strings for gene names.
}
\description{
This internal function is not as accurate (does not retrieve as many gene names as \code{biomaRt}) but can be used without
internet connection.
}
