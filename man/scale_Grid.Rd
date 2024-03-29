% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/miscUtils.R
\name{scale_Grid}
\alias{scale_Grid}
\title{Scale a \code{Grid} object.}
\usage{
scale_Grid(X, fromRange = range_Grid(X), toRange = c(0, 1))
}
\arguments{
\item{X}{An object with class \code{"Grid"} or which can be
coerced into this class.}

\item{fromRange}{A numeric vector of length \code{2} (min and max)
or a matrix with \code{2} rows (min an max) and one column for
each grid dimension in \code{X}. This object gives the original
"old" range.}

\item{toRange}{A numeric vector of length \code{2} (min and max) or
a matrix with \code{2} rows (min an max) and one column for
each grid dimension in \code{X}.This object gives the target "new" range
replacing the old one.}
}
\value{
An object with the same class as \code{X} but with the
each levels rescaled.
}
\description{
Scale a \code{Grid} object.
}
\examples{
myGD <- Grid(levels = list(x = c(1, 3, 10), y = c(0, 1000, 2000, 3000)))
scale_Grid(myGD)
scale_Grid(as.matrix(myGD))
scale_Grid(as.data.frame(myGD))

}
\seealso{
\code{\link{range_Grid}}.
}
