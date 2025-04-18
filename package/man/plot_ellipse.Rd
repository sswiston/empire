% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_ellipse}
\alias{plot_ellipse}
\title{plot_ellipse() function}
\usage{
plot_ellipse(x, y, r, s, a, z = 10)
}
\arguments{
\item{x}{The x coordinate of the ellipse's center point}

\item{y}{The y coordinate of the ellipse's center point}

\item{r}{The r value of the ellipse, which helps define its oblongness (x coordinate of the tilt vector)}

\item{s}{The s value of the ellipse, which helps define its oblongness (y coordinate of the tilt vector)}

\item{a}{The log area of the ellipse}

\item{z}{The height of the tilt vector, default value is 10}
}
\description{
Accepts parameters describing an ellipse and returns a simple plot
}
