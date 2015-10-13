\docType{methods}
\name{aperm}
\alias{aperm}
\title{Permutation of dimensions, or generalized transpose}
\arguments{
  \item{a}{

    Object to be transposed.

  }
  \item{perm}{

    Integer vector giving the permutation.

  }
  \item{...}{

    Additional argument list.

  }
}
\value{
  
  An object with the same class as \code{a} but with
  permuted dimensions.

}
\description{

  Permutation of dimensions, or generalized transpose.

}
\details{
  
  The generalised transposition here is equivalent to the corresponding
  permutation of the columns of the data frame or matrix obtained with
  \code{as.data.frame} or \code{as.matrix}. The order of the rows of the
  data frame or matrix is left unchanged, and is determined by the
  \code{index} slot of the object.  Thus if a vector or matrix of
  response, say \code{Y}, has its elements or rows aligned with the
  design this property remains after the transposition, see
  \bold{Examples}. It must be kept in mind that this "design" order
  can not be guessed when only the factors and their levels are known.

}
\examples{
myGDa <- Grid(nlevels = c("X" = 2, "Y" = 3, "Z" = 4))
myGDb <- aperm(myGDa, perm = c(3, 2, 1))

## evaluation of a function on the permuted array
myGD1 <- Grid(nlevels = c("X" = 2, "Y" = 3, "Z" = 4))
myPerm <- c(3, 2, 1)

myFun <-  function(vec){ sin(vec["X"]) + vec["Y"] - vec["Z"]^2 }

## 'f1' contains the value of the function in the order of 'myGD1'
f1 <- apply_Grid(myGD1, fun = myFun)

## 'f2' contains the value of the function in the order of 'myGD2'
myGD2 <- aperm(myGD1, perm = myPerm)
f2 <- apply_Grid(myGD2, fun = myFun)

## note that 'as.matrix' sorts the observations in the index order
XYZ1 <- as.matrix(myGD1)
XYZ2 <- as.matrix(myGD2)

## check
f = apply(XYZ1, 1, myFun)
cbind(XYZ1, f1 = f1, XYZ2, f2 = f2, f = f)
eps <- sqrt(.Machine$double.eps)
all(abs(f1 - f) < eps) ## should be TRUE
all(abs(f2 - f) < eps) ## should be TRUE    
}
\seealso{
  \code{\link[base:aperm]{aperm}} S3 method from the
  \bold{base} package.
}

