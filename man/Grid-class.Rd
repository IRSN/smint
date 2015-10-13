\name{Grid-class}
\Rdversion{1.1}
\docType{class}
\alias{Grid-class}
\alias{aperm,Grid-method}
\alias{as.data.frame,Grid-method}
\alias{as.matrix,Grid-method}
\alias{dim,Grid-method}
\alias{dimnames,Grid-method}
\alias{dimnames<-,Grid,ANY-method}
\alias{length,Grid-method}
\alias{levels,Grid-method}
\alias{nlevels,Grid-method}
\alias{plot,Grid-method}
\alias{show,Grid-method}

\title{Class \code{"Grid"}}
\description{
 S4 class 
}
\section{Objects from the Class}{
  
  Objects can be created by calls of the form \code{new("Grid",
  ...)}.  or more simply by calling \code{\link{Grid}}.

}
\section{Slots}{
  \describe{
    \item{\code{dim}:}{

      Object of class \code{"integer"} Spatial dimension or number of
      inputs.

    }
    \item{\code{dimNames}:}{

      Object of class \code{"character"} Names of the spatial dimensions
      or inputs.

    }
    \item{\code{levels}:}{

      Object of class \code{"list"} List with length \code{dim}
      containing the levels, i.e. the distinct values taken by
      the corresponding spatial variable.

    }
    \item{\code{index}:}{

      Object of class \code{"array"} Integer index of the combination
      of levels. 

    }
  }
}
\section{Methods}{
  \describe{
    \item{aperm}{

      \code{signature(a = "Grid")}: generalized transposition

    }
    \item{as.data.frame}{

      \code{signature(x = "Grid")}: coercion to a data frame, with
      one row for each combination of the levels.

    }
    \item{as.matrix}{

      \code{signature(x = "Grid")}: coercion to a matrix, with one
      row by combination of the levels.

    }
    \item{checkX}{

      \code{signature(X = "Grid")}: check desing.

    }
    \item{closest}{

      \code{signature(X = "Grid")}: find closest points.

    }
    \item{dim}{

      \code{signature(x = "Grid")}: spatial dimension \eqn{d}.

    }
    \item{dimnames}{

      \code{signature(x = "Grid")}: names of the spatial dimensions.

    }
    \item{dimnames<-}{

      \code{signature(x = "Grid", value = ANY)}: replace the names of the spatial
      dimensions.

    }
    \item{length}{

      \code{signature(x = "Grid")}: number of grid nodes.

    }
    \item{levels}{

      \code{signature(x = "Grid")}: extracts the list of the levels,
      i.e. the possible distinct values for each spatial dimension.

    }
    \item{nlevels}{

      \code{signature(x = "Grid")}: vector of number of levels, one
      value for each spatial dimension.
      
    }
    \item{plot}{
      
      \code{signature(x = "Grid")}: plot an object using
      \code{\link{pairs}}.
      
    }
    \item{sampleIn}{
      
      \code{signature(x = "Grid")}: draw random nodes.
      
    }
    \item{show}{
      
      \code{signature(x = "Grid")}: show grid properties.
      
    }
  }
}
%%\references{}
%%\author{Yves Deville}
%%\note{}
%%\seealso{}
\examples{
showClass("Grid")

## generate a random Grid object, and change its order for tests
set.seed(1234)
GD <- randGrid(nlevels = c("X" = 2, "Y" = 3))
n <- length(GD@index)
ind <- sample(n, size = n, replace = FALSE)
GD@index <- array(ind, dim = nlevels(GD))
X2 <- as.matrix(GD)
GD2 <- as.Grid(X2)
X2Bis <- as.matrix(GD2)
## This hould be zero!
max(abs(X2 - X2Bis))
}
\keyword{classes}
