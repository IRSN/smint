% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gridIntCard.R
\name{gridIntCB}
\alias{gridIntCB}
\title{Grid interpolation in arbitrary dimension through Cardinal
Basis}
\usage{
gridIntCB(
  X,
  Y,
  Xout,
  interpCB = function(x, xout) {
     cardinalBasis_lagrange(x = x, xout = xout)$CB
 },
  intOrder = NULL,
  trace = 1L,
  out_of_bounds = stop,
  ...
)
}
\arguments{
\item{X}{An object that can be coerced into \code{Grid}. This
can be a data.frame or a matrix in \emph{Scattered Data style} in
which case the column number is equal to the spatial dimension
\eqn{d}, and the row number is then equal to the number of nodes
\eqn{n}. But it can also be a \code{Grid} object previously
created.  A data frame or matrix \code{X} will first be coerced
into \code{Grid} by using the the S3 method
\code{\link{as.Grid}}.}

\item{Y}{Response to be interpolated. It must be a vector of
length \eqn{n} equal to the number of nodes. When \code{X} has
class \code{Grid}, the order of the elements in \code{Y} must
conform to the order of the nodes as given in \code{X}, see the
help for \code{\link{gridInt}}.}

\item{Xout}{Interpolation locations. Can be a vector or a
matrix. In the first case, the length of the vector must be equal
to the spatial dimension \eqn{d} as given by \code{xLev} or
\code{X}.  In the second case, each row will be considered as a
response to be interpolated, and the number of columns of
\code{Xout} must be equal to the spatial dimension.}

\item{interpCB}{Function evaluating the interpolation Cardinal
Basis. This function must have as its first 2 formals 'x', and
'xout'.  It must return a matrix with \code{length(x)} columns
and \code{length(xout)} rows. The \eqn{j}-th column is the vector
of the values of the \eqn{j}-th cardinal basis function on
the vector \code{xout}.}

\item{intOrder}{Order of the one-dimensional interpolations. Must
be a permutation of \code{1:d} where \code{d} is the spatial
dimension.  NOT IMPLEMENTED YET. This argument is similar to the
argument of the \code{\link{aperm}} method.}

\item{trace}{Level of verbosity.}

\item{out_of_bounds}{Function to handle Xout outside x (default is stop). Then Xout will be bounded by x range.}

\item{...}{Further arguments to be passed to \code{interpCB}. NOT
IMPLEMENTED YET.}
}
\value{
A single interpolated value if \code{Xout} is a vector or
a row matrix.  If \code{Xout} is a matrix with several rows, the
result is a vector of interpolated values, in the order of the
rows of \code{Xout}.
}
\description{
*****************************************************************************
Grid linear interpolation in arbitrary dimension through Cardinal
Basis.
}
\details{
This is a grid interpolation as in \code{\link{gridInt}} but it is
required here that the one-dimensional interpolation method is
\emph{linear} w.r.t. the vector of interpolated values.  For each
dimension, a one-dimensional interpolation is carried out,
leading to a collection of interpolation problems each with a
dimension reduced by one. The \emph{same cardinal basis} can be
used for all interpolations w.r.t. the same variable, so the
number of call to the \code{interpCB} function is equal to the
interpolation dimension, making this method \emph{very fast} compared to
the general grid interpolation as implemented in
\code{\link{gridInt}}, see the \bold{Examples} section.
}
\examples{
## Natural spline for use through Cardinal Basis in 'gridIntCB'
myInterpCB <-  function(x, xout) cardinalBasis_natSpline(x = x, xout = xout)$CB
## Natural spline for use through Cardinal Basis
myInterp <- function(x, y, xout) {
     spline(x = x, y = y, n = 3 * length(x), method = "natural", xout = xout)$y
}

## generate Grid and function
set.seed(2468)
d <- 5
nLev <- 4L + rpois(d, lambda = 4)
a <- runif(d)
myFun2 <- function(x) exp(-crossprod(a, x^2))
myGD2 <- Grid(nlevels = nLev)
Y2 <- apply_Grid(myGD2, myFun2)
n <- 10
Xout3 <- matrix(runif(n * d), ncol = d)
t1 <- system.time(GI1 <- gridInt(X = myGD2, Y = Y2, Xout = Xout3, interpFun = myInterp))
t2 <- system.time(GI2 <- gridInt(X = myGD2, Y = Y2, Xout = Xout3, interpFun = myInterp,
                                 useC = FALSE))
t3 <- system.time(GI3 <- gridIntCB(X = myGD2, Y = Y2, Xout = Xout3, interpCB = myInterpCB))
df <- data.frame(true = apply(Xout3, 1, myFun2),
                 gridInt_C = GI1, gridInt_R = GI2, gridIntCB = GI3)
head(df)
rbind(gridInt_C = t1, gridInt_C = t2, gridIntCB = t3)

}
\author{
Yves Deville
}
