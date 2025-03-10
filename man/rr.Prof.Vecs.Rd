% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Avg_RelRelief.R
\name{rr.Prof.Vecs}
\alias{rr.Prof.Vecs}
\title{RR of random profiles at all scales}
\usage{
rr.Prof.Vecs(inRas, mn.S, mx.S, stp = 2, samp.P, thresh = 0, sea = NA)
}
\arguments{
\item{inRas}{Input DEM. (matrix NOT SpatRast)}

\item{mn.S}{Min scale.}

\item{mx.S}{Max scale.}

\item{stp}{Scale sequence step. (default = 2)}

\item{samp.P}{Vector of the random profiles sampled.}

\item{thresh}{Threshold elevation to calculate RR above.}

\item{sea}{Direction/orientation of the sea/ocean. ("N", "E", "S", or "W")}
}
\value{
List containing 2 Lists:
Out1 = List of dataframes containing RR for all scales at random profiles
Out2 = List of dataframes containing "moving" Avg of RRs (moving avg. based
scale i.e., avg 3-5, 3-7, 3-9, etc.)
Out 3 = Sample profile.
}
\description{
RR of random profiles at all scales
}
