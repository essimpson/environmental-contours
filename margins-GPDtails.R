library(ismev)
library(evd)

########################################################################
### Function to transform to Laplace margins via ranks and GPD tails ###
########################################################################
laplace.margins <- function(data, p1=0.05, p2=0.95){
  
  ## Apply rank transform
  data.U <- rank(data)/(length(data)+1)
  
  ## Calculate the two thresholds for GPD fitting of the lower/upper tail
  q1 <- quantile(data, p1)
  q2 <- quantile(data, p2)
  
  ## Fit GPDs to both tails
  gpd.fit1 <- gpd.fit(-data, -q1)  
  gpd.fit2 <- gpd.fit(data, q2)
  
  ## Transform to Uniform(0,1)
  data.U[data<q1] <- p1*(1-pgpd(-data[data<q1], loc=-q1, scale=gpd.fit1$mle[1], shape=gpd.fit1$mle[2]))
  data.U[data>q2] <- 1 - (1-p2)*(1-pgpd(data[data>q2], loc=q2, scale=gpd.fit2$mle[1], shape=gpd.fit2$mle[2]))
  
  ## Transform to standard Laplace
  data.L <- data
  data.L[data.U<0.5] <- log(2*data.U[data.U<0.5])
  data.L[data.U>=0.5] <- -log(2*(1-data.U[data.U>=0.5]))
  
  return(data.L)
  
}

##################################################################################################
### Function to back-transform from Laplace margins via ranks with interpolation and GPD tails ###
##################################################################################################
back.transform <- function(x.L, x.orig, p1=0.05, p2=0.95){
  
  ## Deal with the rank-transform first (interpolation to be applied to full distribution and tails subsequently dealt with)
  # First transform from Laplace to Uniform: 
  x.U <- x.L
  x.U[x.L <= 0] <- exp(x.L[x.L<=0])/2
  x.U[x.L > 0]  <- 1 - exp(-x.L[x.L>0])/2
  
  # Original data to Uniform:
  x.orig.U <- rank(x.orig)/(length(x.orig)+1)
  
  # Then interpolate onto original scale:
  x.orig.U.sort <- sort(x.orig.U) 
  x.sort <- sort(x.orig)
  x.new <- NULL
  
  for(i in 1:length(x.U)){
    val <- max(which(x.orig.U.sort<x.U[i]))
    x.new[i] <- x.sort[val] + (x.sort[(val+1)] - x.sort[val])*(x.U[i] - x.orig.U.sort[val])/(x.orig.U.sort[(val+1)] - x.orig.U.sort[val])
  }
  
  ## Fit GPDs to both tails of the original data
  q1 <- quantile(x.orig, p1)
  q2 <- quantile(x.orig, p2)
  gpd.fit1 <- gpd.fit(-x.orig, -q1)  
  gpd.fit2 <- gpd.fit(x.orig, q2)
  
  x.new[x.U<p1] <- -qgpd((1 - x.U[x.U<p1]/p1), loc=-q1, scale=gpd.fit1$mle[1], shape=gpd.fit1$mle[2])
  x.new[x.U>p2] <- qgpd((1 - (1-x.U[x.U>p2])/(1-p2)), loc=q2, scale=gpd.fit2$mle[1], shape=gpd.fit2$mle[2])
  
  return(x.new)
}
