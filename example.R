rm(list=ls())

### Load required functions from other files:
source("estimation-Contour1.R") # Functions for estimating environmental contours for definition 1
source("estimation-Contour2.R") # Functions for estimating environmental contours for definition 2
source("margins-GPDtails.R")    # Functions for transforming to standard Laplace margins

### Load required packages:
library(scales)
library(sp)
library(VGAM)

#################################################################

# Example - Independent Laplace variables
data <- matrix(rlaplace(365*30*2), ncol=2)

# Estimating a contour where the proportion of mass outside is 1%
# Contour definition 1
contour1 <- contourEstimation1(data, quant=0.99, thresh.quant=0.5, nei=100, len=200, knot.num=24)

plot(data, 
     main="Independence: Contour 1", xlab=expression(X[L]), ylab=expression(Y[L]), 
     cex.lab=1.5, cex.main=1.5, cex.axis=1.5, 
     asp=1, col=alpha(1,0.3), pch=16, xlim=c(-12,12), ylim=c(-12,12))
points(contour1$contour, col=6, type="l", lwd=3)


# Contour definition 2
contour2 <- contourEstimation2(data, quant=0.99, thresh.quant=0.5, nei=100, len=200, knot.num=24)

plot(data, 
     main="Independence: Contour 2", xlab=expression(X[L]), ylab=expression(Y[L]), 
     cex.lab=1.5, cex.main=1.5, cex.axis=1.5, 
     asp=1, col=alpha(1,0.3), pch=16, xlim=c(-12,12), ylim=c(-12,12))
points(contour2$contour, col=6, type="l", lwd=3, lty=2)
