################################################################################
# Updated version of the R code for the analysis in:
#
#   "A penalized framework for distributed lag non-linear models"
#   Biometrics, 2017
#   Antonio Gasparrini, Fabian Scheipl, Ben Amstrong, and Michael G. Kenward
#   http://www.ag-myresearch.com/2017_gasparrini_biomet.html
#
# Update: 15 Jan 2017
# * an up-to-date version of this code is available at:
#   https://github.com/gasparrini/2017_gasparrini_biomet_Rcodedata
################################################################################

################################################################################
# FIRST EXAMPLE
# LOAD THE DATA AND RUN THE MODELS (EXTERNAL METHOD)
################################################################################

# LOAD PACKAGES AND FUNCTIONS
library(dlnm) ; library(mgcv) ; library(splines) ; library(tsModel)

# LOAD DATA
london <- read.csv("london.csv")

################################################################################
# GLM WITH KNOTS SPECIFIED A PRIORI (GASPARRINI BMCmrm 2014)

# DEFINE THE CROSS-BASIS
vk <- equalknots(london$tmean,nk=2)
lk <- logknots(25,nk=3)
cbglm1 <- crossbasis(london$tmean, lag=25, argvar=list(fun="bs",degree=2,
  knots=vk), arglag=list(knots=lk))
summary(cbglm1)

# RUN THE MODEL AND PREDICT
library(splines)
glm1 <- glm(death~cbglm1+ns(time,10*14)+dow,family=quasipoisson(),london)
pred3dglm1 <- crosspred(cbglm1,glm1,at=-3:29,cen=20)
predslglm1 <- crosspred(cbglm1,glm1,by=0.2,bylag=0.2,cen=20)

# PLOTS
plot(pred3dglm1,xlab="Temperature (C)",zlab="RR",zlim=c(0.88,1.45),xlim=c(-5,30),
  ltheta=170,phi=35,lphi=30,main="Original GLM")
plot(predslglm1,"overall",ylab="RR",xlab="Temperature (C)",xlim=c(-5,30),
  ylim=c(0.5,3.5),lwd=1.5,main="Original GLM")
plot(predslglm1,var=29,xlab="Lag (days)",ylab="RR",ylim=c(0.9,1.4),lwd=1.5,
  main="Original GLM")

################################################################################
# GLM WITH AIC-BASED KNOT SELECTION

# Q-AIC FUNCTION
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}

# KNOTS GRID
grid <- as.matrix(expand.grid(var=1:8,lag=1:8))

# SEARCH (TAKES ~40sec IN A 2.4 GHz PC)
system.time({
aicval <- sapply(seq(nrow(grid)), function(i) {
  vk <- equalknots(london$tmean,nk=grid[i,1])
  lk <- logknots(25,nk=grid[i,2])
  cb <- crossbasis(london$tmean,lag=25,argvar=list(fun="bs",degree=2,knots=vk),
    arglag=list(knots=lk))
  m <- glm(death~cb+ns(time,10*14)+dow,family=quasipoisson(),london)
  return(fqaic(m))
})
})

# BEST FITTING MODEL
(best <- grid[which.min(aicval),])
plot(aicval,col=2,pch=19)

# DEFINE THE CROSS-BASIS
vk <- equalknots(london$tmean,nk=best[1])
lk <- logknots(25,nk=best[2])
cbglm2 <- crossbasis(london$tmean,lag=25,argvar=list(fun="bs",degree=2,
  knots=vk),arglag=list(knots=lk))
summary(cbglm2)

# RUN THE MODEL AND PREDICT
glm2 <- glm(death~cbglm2+ns(time,10*14)+dow,family=quasipoisson(),london)
pred3dglm2 <- crosspred(cbglm2,glm2,at=-3:29,cen=20)
predslglm2 <- crosspred(cbglm2,glm2,by=0.2,bylag=0.2,cen=20)

# PLOTS
plot(pred3dglm2,xlab="Temperature (C)",zlab="RR",zlim=c(0.88,1.45),xlim=c(-5,30),
  ltheta=170,phi=35,lphi=30,main="GLM with AIC-based knot selection")
plot(predslglm2,"overall",ylab="RR",xlab="Temperature (C)",xlim=c(-5,30),
  ylim=c(0.5,3.5),lwd=1.5,main="GLM with AIC-based knot selection")
plot(predslglm2,var=29,xlab="Lag (days)",ylab="RR",ylim=c(0.9,1.4),lwd=1.5,
  main="GLM with AIC-based knot selection")

################################################################################
# GAM WITH DEFAULT PENALTIES

# DEFINE THE CROSS-BASIS
# NB: df IN argvar SET TO 9, AS INTERCEPT IS EXCLUDED AUTOMATICALLY
# (FOR COMPATIBILITY WITH INTERNAL METHOD)
cbgam1 <- crossbasis(london$tmean,lag=25,argvar=list(fun="ps",df=9),
  arglag=list(fun="ps",df=10))
summary(cbgam1)

# DEFINE THE PENALTY MATRICES
cbgam1Pen <- cbPen(cbgam1)

# RUN THE GAM MODEL AND PREDICT (TAKES ~34sec IN A 2.4 GHz PC)
system.time({
gam1 <- gam(death~cbgam1+ns(time,10*14)+dow,family=quasipoisson(),london,
  paraPen=list(cbgam1=cbgam1Pen), method='REML')
})
pred3dgam1 <- crosspred(cbgam1,gam1,at=-3:29,cen=20)
predslgam1 <- crosspred(cbgam1,gam1,by=0.2,bylag=0.2,cen=20)

# CHECK CONVERGENCE, SMOOTHING PARAMETERS AND EDF
gam1$converged
gam1$sp
sum(gam1$edf[2:91])

# PLOTS
plot(pred3dgam1,xlab="Temperature (C)",zlab="RR",zlim=c(0.88,1.45),xlim=c(-5,30),
 ltheta=170,phi=35,lphi=30,main="GAM with default penalties")
plot(predslgam1,"overall",ylab="RR",xlab="Temperature (C)",xlim=c(-5,30),
  ylim=c(0.5,3.5),lwd=1.5,main="GAM with default penalties")
plot(predslgam1,var=29,xlab="Lag (days)",ylab="RR",ylim=c(0.9,1.4),lwd=1.5,
  main="GAM with default penalties")

################################################################################
# GAM WITH DOUBLY VARYING PENALTY ON THE LAG

# DEFINE THE CROSS-BASIS
# PS: EXCLUDE THE DEFAULT PENALTY FROM THE LAG-RESPONSE FUNCTION WITH fx=T
cbgam2 <- crossbasis(london$tmean, lag=25, argvar=list(fun="ps",df=9),
  arglag=list(fun="ps",df=10,fx=T))
summary(cbgam2)

# DEFINE THE DOUBLY VARYING PENALTY MATRICES
# VARYING DIFFERENCE PENALTY APPLIED TO LAGS (EQ. 8b)
C <- do.call('onebasis',c(list(x=0:25,fun="ps",df=10,intercept=T)))
D <- diff(diag(25+1),diff=2)
P <- diag((seq(0,25-2))^2)
Slag1 <- t(C) %*% t(D) %*% P %*% D %*% C
# VARYING RIDGE PENALTY APPLIED TO COEFFICIENTS (Eq. 7a)
Slag2 <- diag(rep(0:1,c(6,4)))
cbgam2Pen <- cbPen(cbgam2,addSlag=list(Slag1,Slag2))

# RUN THE GAM MODEL AND PREDICT (TAKES ~14sec IN A 2.4 GHz PC)
system.time({
gam2 <- gam(death~cbgam2+ns(time, 10*14)+dow,family=quasipoisson(),london,
  paraPen=list(cbgam2=cbgam2Pen), method='REML')
})
pred3dgam2 <- crosspred(cbgam2,gam2,at=-3:29,cen=20)
predslgam2 <- crosspred(cbgam2,gam2,by=0.2,bylag=0.2,cen=20)

# CHECK CONVERGENCE, SMOOTHING PARAMETERS AND EDF
gam2$converged
gam2$sp
sum(gam2$edf[2:91])

# PLOTS
plot(pred3dgam2,xlab="Temperature (C)",zlab="RR",zlim=c(0.88,1.45),xlim=c(-5,30),
 ltheta=170,phi=35,lphi=30,main="GAM with doubly varying penalties")
plot(predslgam2,"overall",ylab="RR",xlab="Temperature (C)",xlim=c(-5,30),
  ylim=c(0.5,3.5),lwd=1.5,main="GAM with doubly varying penalties")
plot(predslgam2,var=29,xlab="Lag (days)",ylab="RR",ylim=c(0.9,1.4),lwd=1.5,
  main="GAM with doubly varying penalties")

################################################################################
# GAM WITH ONE DIMENSION UNPENALIZED BY USING PARAMETRIC FUNCTIONS

# DEFINE THE CROSS-BASIS
# NB: DOUBLE THRESHOLD AS EXPOSURE-RESPONSE
cbgam3 <- crossbasis(london$tmean, lag=25, argvar=list(fun="thr",thr=c(17,21)),
  arglag=list(fun="ps",df=10))
summary(cbgam3)

# USE THE PENALTY MATRICES DEFINE ABOVE
cbgam3Pen <- cbPen(cbgam3,addSlag=list(Slag1,Slag2))

# RUN THE GAM MODEL AND PREDICT (TAKES ~75sec IN A 2.4 GHz PC)
system.time({
gam3 <- gam(death~cbgam3+ns(time, 10*14)+dow,family=quasipoisson(),london,
  paraPen=list(cbgam3=cbgam3Pen), method='REML')
})
pred3dgam3 <- crosspred(cbgam3,gam3,at=-3:29)
predslgam3 <- crosspred(cbgam3,gam3,by=0.2,bylag=0.2)

# CHECK CONVERGENCE, SMOOTHING PARAMETERS AND EDF
gam3$converged
gam3$sp
sum(gam3$edf[2:21])

# PLOTS
plot(pred3dgam3,xlab="Temperature (C)",zlab="RR",zlim=c(0.88,1.45),xlim=c(-5,30),
  ltheta=170,phi=35,lphi=30,main="Mix of penalized and unpenalized")
plot(predslgam3,"overall",ylab="RR",xlab="Temperature (C)",xlim=c(-5,30),
  ylim=c(0.8,2.2),lwd=1.5,main="Mix of penalized and unpenalized")
plot(predslgam3,var=29,xlab="Lag (days)",ylab="RR",ylim=c(0.9,1.4),lwd=1.5,
  main="Mix of penalized and unpenalized")

#
