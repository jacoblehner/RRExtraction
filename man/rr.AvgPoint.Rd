% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Avg_RelRelief.R
\name{rr.AvgPoint}
\alias{rr.AvgPoint}
\title{Multi-scale average RR at a given cell}
\usage{
rr.AvgPoint(inRas, inMn = 7, inMx = 21, stp = 2, ii, jj)
}
\arguments{
\item{inRas}{Input DEM. (SpatRast)}

\item{inMn}{Minimum scale. (default = 7)}

\item{inMx}{Maximum scale. (default = 21)}

\item{stp}{Scale step. (default = 2)}

\item{ii}{Row position of cell.}

\item{jj}{Column position of cell.}
}
\value{
Average relative relief at chosen cell. (numeric)
}
\description{
Multi-scale average RR at a given cell
}
