% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ellipse.R
\name{get_ellipse}
\alias{get_ellipse}
\title{get_ellipse() function}
\usage{
get_ellipse(data = NULL, z = 10, confidence = 0.95)
}
\arguments{
\item{data}{The spatial POLYGON or MULTIPOLYGON to turn into an ellipse. Alternatively, the occurrence points to turn into an ellipse, as a data.frame or MULTIPOINT object.}

\item{z}{The height of the tilt vector for generating the ellipses, default value is 10}

\item{confidence}{The level of confidence used to generate the confidence ellipse, default value is 0.95}
}
\description{
Takes a spatial POLYGON or MULTIPOLYGON (sf package) and generates ellipse parameters x, y, r, s, and a. Returns a vector (x,y,r,s,a,z). This is done by sampling points from within the spatial polygon/multipolygon, then creating a confidence ellipse around those points using the stats function cov.wt(). Alternatively, can provide a set of occurrence points in a data.frame or MULTIPOINT (sf) object.
}
