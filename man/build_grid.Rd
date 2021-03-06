% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algorithms.R
\name{build_grid}
\alias{build_grid}
\title{Build the grid for the grid search algorithm with a zoom.}
\usage{
build_grid(...)
}
\arguments{
\item{...}{a sequence of vectors or lists containing the information about the grid to be built, see Usage and Details.}
}
\value{
a new object of the class GRID with the grid ready for the grid search with a zoom.

The object contains the following components:
\item{grid}{the grid}
\item{size}{number of points in the grid}
\item{npar}{number of arguments or parameters}
}
\description{
This function builds the grid for the grid search algorithm with a zoom.
}
\details{
The argument \code{...} is a sequence of vectors or lists contatining the information about the grid to be built.
Each element in the sequence is either a vector or a list taking one of the following forms

- \code{x}, if \code{x} is already an sequence of the grid points for the corresponding argument.

- \code{c(from=, to=, by=)}

- \code{c(from=, to=, length=)}

- \code{list(from=, to=, by=)}

- \code{list(from=, to=, length=)}

where

- \code{from}: the min of the argument of the target function

- \code{to}: the max of the argument of the target function

- \code{by}: the increment of the sequence

- \code{length}: desired length.

There are many different ways to organize the points on the grid for certain argument of the target function,
the user can make them freely and input directly by \code{build_grid(x, ...)}.
Notice that \code{x} does not need to be increasing, as the function will sort it.
The design that \code{x} does not need to be increasing makes it convenient for the user
to interpolate more points at some region without considering to sort it all the time.

When \code{by} is provided, the \code{length} will be ignored.
So if the user wanna specify the \code{length}, please do not use \code{by}.

The order of the sequence \code{...} matters as it represents the order of the corresponding arguments of the target function to be optimized.
}
\examples{
vx = 1:5
build_grid(vx, c(from=1, to=2, by=.2), list(from=3, to=4, length=5))

}
\seealso{
\code{\link{grid_search_check}}, \code{\link{grid_search}}
}
\author{
Yukai Yang, \email{yukai.yang@statistik.uu.se}
}
\keyword{algorithms}
