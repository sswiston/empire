% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_prob.R
\name{get_expectations}
\alias{get_expectations}
\title{get_expectations() function}
\usage{
get_expectations(tree, param, XY, W, V)
}
\arguments{
\item{dataTree}{A dataTree object containing a detailedTree object and dataframes V and W containing values for cladogenetic scenarios at nodes}
}
\description{
Accepts a dataTree object and returns the expected values of x and y values at the tips based on the root values, tree structure, and data-augmented values
}
