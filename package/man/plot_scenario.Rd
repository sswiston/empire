% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_scenario}
\alias{plot_scenario}
\title{plot_scenario() function}
\usage{
plot_scenario(x, y, r, s, a, d, m, c, h, z = 10, alpha = -3)
}
\arguments{
\item{x}{The x coordinate of the ancestral ellipse's center point}

\item{y}{The y coordinate of the ancestral ellipse's center point}

\item{r}{The r value of the ancestral ellipse, which helps define its oblongness (x coordinate of the tilt vector)}

\item{s}{The s value of the ancestral ellipse, which helps define its oblongness (y coordinate of the tilt vector)}

\item{a}{The log area of the ancestral ellipse}

\item{d}{The daughter configuration during cladogenesis (0,1), where 0 makes the left daughter D1 and 1 makes the right daughter D1}

\item{m}{The cladogenetic mode (0,1), where 0 is sympatric and 1 is allopatric}

\item{c}{The c value (not index) of the cladogenetic scenario's concentric circle --> ex. 0.5}

\item{h}{The h value (not index) of the cladogenetic scenario's direction line --> ex. pi/4, if required}

\item{z}{The height of the tilt vector, default value is 10}

\item{alpha}{The alpha parameter, default value is -3}
}
\description{
Accepts parameters describing an ellipse and cladogenetic scenario and returns a simple plot
}
