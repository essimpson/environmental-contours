library(evgam)
library(mgcv)

source("margins-GPDtails.R")

###########################################################################################################
### Function to carry out boundary and parameter estimation with the spline degree automatically chosen ###
###########################################################################################################
contourEstimation1 <- function(data, quant=0.999, thresh.quant=0.50, nei=100, len=200, knot.num=24){
  
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
  
  # For each value of W*: calculate the largest angular distance (from W*) in each neighbourhood; 
  # calculate estimates of the radial quantiles from a standard GPD fit with empirical quantile threshold:
  thresh.move <- NULL
  R.W.quant.ind <- NULL
  
  for(sec.num in 1:length(Wstar)){
    sec.W <- Wstar[sec.num]
    eps.W <- sort(abs(W.L - sec.W))[nei]
    nei.W <- which(abs(W.L - sec.W) <= eps.W)
    
    thresh.move[sec.num] <- quantile(R.L[nei.W],thresh.quant)
    fit <- ismev::gpd.fit(R.L[nei.W], threshold=thresh.move[sec.num], npy=length(nei.W))
    R.W.quant.ind[sec.num] <- evd::qgpd((quant-thresh.quant)/(1-thresh.quant), loc=thresh.move[sec.num], scale=fit$mle[1], shape=fit$mle[2])
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
  R.matrix <- as.data.frame(cbind(W.L[exceeds1==T], (R.L[exceeds1==T]-thresh.smooth1[exceeds1==T])))
  names(R.matrix) <- c("W.L","R.L") 
  spl=paste('R.L ~ s(W.L, bs="cp", k=',knot.num,', m=c(0,2))')
  fmla_gpd <- list(as.formula(spl), ~1)
  m_gpd <- evgam(fmla_gpd, data=R.matrix, family="gpd", knots=knots)
  gpd_pred <- predict(m_gpd, newdata=R.matrix.pred, prob=(quant-thresh.quant)/(1-thresh.quant))
  R.W.quant1 <- thresh.smooth1a + gpd_pred$`q`
  errors <- sum(abs(R.W.quant1 - R.W.quant.ind))
  
  R.matrix <- as.data.frame(cbind(W.L[exceeds2==T], (R.L[exceeds2==T]-thresh.smooth2[exceeds2==T])))
  names(R.matrix) <- c("W.L","R.L") 
  spl=paste('R.L ~ s(W.L, bs="cp", k=',knot.num,', m=c(1,2))')
  fmla_gpd <- list(as.formula(spl), ~1)
  m_gpd <- evgam(fmla_gpd, data=R.matrix, family="gpd", knots=knots)
  gpd_pred <- predict(m_gpd, newdata=R.matrix.pred, prob=(quant-thresh.quant)/(1-thresh.quant))
  R.W.quant2 <- thresh.smooth2a + gpd_pred$`q`
  errors[2] <- sum(abs(R.W.quant2 - R.W.quant.ind))
  
  R.matrix <- as.data.frame(cbind(W.L[exceeds3==T], (R.L[exceeds3==T]-thresh.smooth3[exceeds3==T])))
  names(R.matrix) <- c("W.L","R.L") 
  spl=paste('R.L ~ s(W.L, bs="cp", k=',knot.num,', m=c(2,2))')
  fmla_gpd <- list(as.formula(spl), ~1)
  m_gpd <- evgam(fmla_gpd, data=R.matrix, family="gpd", knots=knots)
  gpd_pred <- predict(m_gpd, newdata=R.matrix.pred, prob=(quant-thresh.quant)/(1-thresh.quant))
  R.W.quant3 <- thresh.smooth3a + gpd_pred$`q`
  errors[3] <- sum(abs(R.W.quant3 - R.W.quant.ind))
  
  # Select R-quantile estimates using spline degree with smallest error  
  degree <- which.min(errors)
  R.W.quant <- R.W.quant1*(degree==1) + R.W.quant2*(degree==2) + R.W.quant3*(degree==3)

  # Transform back to original coordinates to obtain contour estimates:
  contour.L <- cbind(R.W.quant*cos(Wstar),R.W.quant*sin(Wstar))
  contour.orig <- cbind(back.transform(contour.L[,1], data[,1]), back.transform(contour.L[,2], data[,2]))

  return(list(data.L=data.L, contour=contour.orig, contour.Laplace=contour.L, degree=degree, R.L=R.L, W.L=W.L))
  
}


