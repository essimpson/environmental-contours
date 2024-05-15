library(evgam)
library(mgcv)

source("margins-GPDtails.R")

####################################################################
### Function to integrate over for scaling constant in quantiles ###
####################################################################
scaling.fun <- function(w, eps=eps, W.L.dens){
  
  f.W <- NULL
  
  for(i in 1:length(w)){
    val <- max(which(W.L.dens$x<w[i]))
    f.W[i] <- W.L.dens$y[val] + (W.L.dens$y[(val+1)] - W.L.dens$y[val])*(w[i] - W.L.dens$x[val])/(W.L.dens$x[(val+1)] - W.L.dens$x[val])
  }
  f.W <- 3*f.W
  
  return(apply(cbind(rep(1,length(w)),f.W/eps), 1, min))

}
  
###############################################################
### Function to carry out boundary and parameter estimation ###
###############################################################
contourEstimation2 <- function(data, quant=0.999, thresh.quant=0.50, nei=100, len=200, knot.num=24){
  
  eps <- (1-quant)/(2*pi)
  
  # Fix knot locations
  if(knot.num%%8!=0){return("The specified number of knots should be a multiple of 8.")}
  knots <- list(W.L=seq(-pi, pi, length.out=knot.num+1))
  
  # Transform to Laplace margins:
  data.L <- apply(data, 2, laplace.margins)
  
  # Calculate polar coordinates:
  R.L <- sqrt(data.L[,1]^2  + data.L[,2]^2)
  W.L <- atan2(data.L[,2],data.L[,1])

  # Take log(R) for threshold calculations:
  logR.L <- log(R.L)
  
  # Vector of angles at which to estimate the boundary:
  if(len%%8!=0){return("The length of Wstar should be a multiple of 8.")}
  Wstar <- seq(-pi,pi,length.out=(len+1))
  
  # Kernel density estimation for the angles evaluated at W* and original W values
  # Repeat data to ensure cyclic consistency - multiply densities by 3 in required range
  W.L.rep <- c((W.L-2*pi),W.L,(W.L+2*pi))
  W.L.dens <- density(W.L.rep, bw="SJ")
  f.Wstar <- NULL
  for(i in 1:length(Wstar)){
    val <- max(which(W.L.dens$x<Wstar[i]))
    f.Wstar[i] <- W.L.dens$y[val] + (W.L.dens$y[(val+1)] - W.L.dens$y[val])*(Wstar[i] - W.L.dens$x[val])/(W.L.dens$x[(val+1)] - W.L.dens$x[val])
  }
  f.Wstar <- 3*f.Wstar
  
  # Radial quantiles of interest
  int <- integrate(scaling.fun, -pi, pi, eps=eps, W.L.dens=W.L.dens)$value
  con <- (1-quant)/int
  q.Wstar <- 1 - (con/apply(cbind(rep(eps,length(f.Wstar)),f.Wstar),1,max))

  # For each value of W*: calculate the largest angular distance (from W*) in each neighbourhood; 
  # calculate estimates of the radial quantiles from a standard GPD fit with empirical quantile threshold:
  thresh.move <- NULL
  R.W.quant.ind <- rep(NA, length(Wstar))
  
  for(sec.num in 1:length(Wstar)){
    sec.W <- Wstar[sec.num]
    eps.W <- sort(abs(W.L - sec.W))[nei]
    nei.W <- which(abs(W.L - sec.W) <= eps.W)
    
    thresh.move[sec.num] <- quantile(R.L[nei.W],thresh.quant)
    if(q.Wstar[sec.num]<=thresh.quant){
      R.W.quant.ind[sec.num] <- quantile(R.L[nei.W],q.Wstar[sec.num])
    }
    if(q.Wstar[sec.num]>thresh.quant){
      fit <- ismev::gpd.fit(R.L[nei.W], threshold=thresh.move[sec.num], npy=length(nei.W))
      R.W.quant.ind[sec.num] <- evd::qgpd((q.Wstar[sec.num]-thresh.quant)/(1-thresh.quant), loc=thresh.move[sec.num], scale=fit$mle[1], shape=fit$mle[2])
    }
  }
  
  # Use asymmetric Laplace for quantile regression of threshold:
  # Fit on log(R) then back-transform:
  fmla_ald <- paste('logR.L ~ s(W.L, bs="cp", k=',knot.num,', m=c(0,2))')
  m_ald <- evgam(as.formula(fmla_ald), data=as.data.frame(cbind(logR.L,W.L)), family="ald", ald.args=list(tau=thresh.quant), knots=knots)
  thresh.smooth1 <- exp(predict(m_ald, newdata=list(W.L=W.L))$location)
  thresh.smooth1a <- exp(predict(m_ald, newdata=list(W.L=Wstar))$location)
  
  fmla_ald <- paste('logR.L ~ s(W.L, bs="cp", k=',knot.num,', m=c(1,2))')
  m_ald <- evgam(as.formula(fmla_ald), data=as.data.frame(cbind(logR.L,W.L)), family="ald", ald.args=list(tau=thresh.quant), knots=knots)
  thresh.smooth2 <- exp(predict(m_ald, newdata=list(W.L=W.L))$location)
  thresh.smooth2a <- exp(predict(m_ald, newdata=list(W.L=Wstar))$location)
  
  fmla_ald <- paste('logR.L ~ s(W.L, bs="cp", k=',knot.num,', m=c(2,2))')
  m_ald <- evgam(as.formula(fmla_ald), data=as.data.frame(cbind(logR.L,W.L)), family="ald", ald.args=list(tau=thresh.quant), knots=knots)
  thresh.smooth3 <- exp(predict(m_ald, newdata=list(W.L=W.L))$location)
  thresh.smooth3a <- exp(predict(m_ald, newdata=list(W.L=Wstar))$location)
  
  # Extract threshold exceedances for subsequent GPD fit:
  exceeds1 <- apply(cbind(R.L,thresh.smooth1),1,function(x){x[1]>x[2]})
  exceeds2 <- apply(cbind(R.L,thresh.smooth2),1,function(x){x[1]>x[2]})
  exceeds3 <- apply(cbind(R.L,thresh.smooth3),1,function(x){x[1]>x[2]})
  
  R.matrix.pred <- as.data.frame(cbind(Wstar, rep(NA,length(Wstar)))) 
  names(R.matrix.pred) <- c("W.L","R.L") 

  # Fit the GPD-GAM (with a spline in the scale parameter and constant shape parameter)
  # to the radial components, using the W* values as a covariate in the spline.
  # Calculate the required radial quantile for each W*.
  # Do this for splines of degree 1,2,3 and select the best by comparing to the individual GPD estimates:
  R.W.quant1 <- rep(NA, length(Wstar))
  for(iter in 1:length(Wstar)){
    if(q.Wstar[iter]<=thresh.quant){
      fmla_ald <- paste('logR.L ~ s(W.L, bs="cp", k=',knot.num,', m=c(0,2))')
      m_ald <- evgam(as.formula(fmla_ald), data=as.data.frame(cbind(logR.L,W.L)), family="ald", ald.args=list(tau=q.Wstar[iter]), knots=knots)
      R.W.quant1[iter] <- exp(predict(m_ald, newdata=list(W.L=Wstar))$location)[iter]
    }
  }
  R.matrix <- as.data.frame(cbind(W.L[exceeds1==T], (R.L[exceeds1==T]-thresh.smooth1[exceeds1==T])))
  names(R.matrix) <- c("W.L","R.L") 
  spl=paste('R.L ~ s(W.L, bs="cp", k=',knot.num,', m=c(0,2))')
  fmla_gpd <- list(as.formula(spl), ~1)
  m_gpd <- evgam(fmla_gpd, data=R.matrix, family="gpd", knots=knots)
  gpd_pred <- predict(m_gpd, newdata=R.matrix.pred, prob=(q.Wstar-thresh.quant)/(1-thresh.quant))
  R.W.quant1[q.Wstar>thresh.quant] <- thresh.smooth1a[q.Wstar>thresh.quant] + diag(as.matrix(gpd_pred))[q.Wstar>thresh.quant]
  errors <- sum(abs(R.W.quant1 - R.W.quant.ind))

  R.W.quant2 <- rep(NA, length(Wstar))
  for(iter in 1:length(Wstar)){
    if(q.Wstar[iter]<=thresh.quant){
      fmla_ald <- paste('logR.L ~ s(W.L, bs="cp", k=',knot.num,', m=c(1,2))')
      m_ald <- evgam(as.formula(fmla_ald), data=as.data.frame(cbind(logR.L,W.L)), family="ald", ald.args=list(tau=q.Wstar[iter]), knots=knots)
      R.W.quant2[iter] <- exp(predict(m_ald, newdata=list(W.L=Wstar))$location)[iter]
    }
  }
  R.matrix <- as.data.frame(cbind(W.L[exceeds2==T], (R.L[exceeds2==T]-thresh.smooth2[exceeds2==T])))
  names(R.matrix) <- c("W.L","R.L") 
  spl=paste('R.L ~ s(W.L, bs="cp", k=',knot.num,', m=c(1,2))')
  fmla_gpd <- list(as.formula(spl), ~1)
  m_gpd <- evgam(fmla_gpd, data=R.matrix, family="gpd", knots=knots)
  gpd_pred <- predict(m_gpd, newdata=R.matrix.pred, prob=(q.Wstar-thresh.quant)/(1-thresh.quant))
  R.W.quant2[q.Wstar>thresh.quant] <- thresh.smooth2a[q.Wstar>thresh.quant] + diag(as.matrix(gpd_pred))[q.Wstar>thresh.quant]
  errors[2] <- sum(abs(R.W.quant2 - R.W.quant.ind))

  R.W.quant3 <- rep(NA, length(Wstar))
  for(iter in 1:length(Wstar)){
    if(q.Wstar[iter]<=thresh.quant){
      fmla_ald <- paste('logR.L ~ s(W.L, bs="cp", k=',knot.num,', m=c(2,2))')
      m_ald <- evgam(as.formula(fmla_ald), data=as.data.frame(cbind(logR.L,W.L)), family="ald", ald.args=list(tau=q.Wstar[iter]), knots=knots)
      R.W.quant3[iter] <- exp(predict(m_ald, newdata=list(W.L=Wstar))$location)[iter]
    }
  }
  R.matrix <- as.data.frame(cbind(W.L[exceeds3==T], (R.L[exceeds3==T]-thresh.smooth3[exceeds3==T])))
  names(R.matrix) <- c("W.L","R.L") 
  spl=paste('R.L ~ s(W.L, bs="cp", k=',knot.num,', m=c(2,2))')
  fmla_gpd <- list(as.formula(spl), ~1)
  m_gpd <- evgam(fmla_gpd, data=R.matrix, family="gpd", knots=knots)
  gpd_pred <- predict(m_gpd, newdata=R.matrix.pred, prob=(q.Wstar-thresh.quant)/(1-thresh.quant))
  R.W.quant3[q.Wstar>thresh.quant] <- thresh.smooth3a[q.Wstar>thresh.quant] + diag(as.matrix(gpd_pred))[q.Wstar>thresh.quant]
  errors[3] <- sum(abs(R.W.quant3 - R.W.quant.ind))
  
  # Select R-quantile estimates using spline degree with smallest error  
  degree <- which.min(errors)
  R.W.quant <- R.W.quant1*(degree==1) + R.W.quant2*(degree==2) + R.W.quant3*(degree==3)

  # Transform back to original coordinates to obtain boundary estimates:
  contour.L <- cbind(R.W.quant*cos(Wstar),R.W.quant*sin(Wstar))
  contour.orig <- cbind(back.transform(contour.L[,1], data[,1]), back.transform(contour.L[,2], data[,2]))

  return(list(data.L=data.L, contour=contour.orig, contour.Laplace=contour.L, degree=degree, R.L=R.L, W.L=W.L))
  
}
