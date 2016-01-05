# smint: Smooth Multivariate Interpolation for Gridded and Scattered Data

Installation
------------

You can install the latest version of the code using the `devtools` R package.

```
# Install devtools, if you haven't already.
install.packages("devtools")

library(devtools)
install_github("IRSN/smint")
```

Usage 
-----------
Complete user guide is here: [smintGuide.pdf](https://github.com/IRSN/smint/blob/ec795cf4457437972285bd435d546f62427cbeff/inst/doc/smintGuide.pdf)

Basic example:
```r
     set.seed(12345)
     n <- 6; nout <- 300L
     x <- sort(runif(n))
     xout <- sort(runif(nout, min = x[1], max = x[n]))
     y <- sin(2 * pi * x)
     cI0 <- interp_ceschino(x = x, xout = xout, y = y)
     
     ## compare with a natural spline
     require(splines)
     spI <- interpSpline(x, y)
     spPred0 <- predict(spI, xout)
     
     plot(xout, sin(2 * pi * xout), type = "l", col = "black", lwd = 2,
          xlab = "x", ylab = "f(x)", main = "Interpolations")
     abline(v = x, col = "gray")
     lines(cI0, type = "l", col = "SpringGreen3", lty = 2, lwd = 2)
     lines(spPred0, type = "l", col = "SteelBlue2", lty = 3, lwd = 2)
     points(x, y, type = "p", pch = 21, col = "red",
            bg = "yellow", lwd = 2)
     legend("topright", legend = c("true", "Ceschino", "nat. spline"),
             col = c("black", "SpringGreen3", "SteelBlue2"),
             lty = 1:3, lwd = rep(2, 3))
```
![](Rplot.png)
