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
# REPRODUCE THE GRAPHS
################################################################################

# PARAMETERS FOR COLOURED GRAPHS
border <- NULL
col3d <- "lightskyblue"
col <- 2

# PARAMETERS FOR BLACK&WHITE GRAPHS (UNCOMMENT IF NEEDED)
# border <- grey(0.3)
# col3d <- grey(0.99)
# col <- 1

################################################################################

pdf("fig3.pdf",height=12,width=9)
layout(matrix(1:12,ncol=3,byrow=TRUE))

################################################################################
# MAIN GAM (DOUBLY VARYING PENALTY ON THE LAG)

par(mar=c(1,1,3,1))

plot(pred3dgam2,xlab="Temperature (C)",zlab="RR",zlim=c(0.9,1.4),xlim=c(-5,30),
  ltheta=170,phi=35,lphi=30,shade=0.75,cex.axis=0.7,cex.lab=1,
  border=border,col=col3d,main="Exposure-lag-response")
mtext("GAM with doubly varying penalty",cex=0.6)

par(mar=c(5,4,4,1))

plot(predslgam2,"overall",ylab="RR",xlab="Temperature (C)",xlim=c(-5,30),
  ylim=c(0.5,3.5),lwd=1.5,col=col,main="Overall cumulative exposure-response")
mtext("GAM with doubly varying penalty",cex=0.6)

plot(predslgam2,var=29,xlab="Lag (days)",ylab="RR",ylim=c(0.9,1.4),
  lwd=1.5,col=col,main="Lag-response at 29C")
mtext("GAM with doubly varying penalty",cex=0.6)

################################################################################
# GLM WITH KNOTS SPECIFIED A PRIORI (GASPARRINI BMCmrm 2014)

par(mar=c(1,1,3,1))

plot(pred3dglm1,xlab="Temperature (C)",zlab="RR",zlim=c(0.9,1.4),xlim=c(-5,30),
  ltheta=170,phi=35,lphi=30,shade=0.75,cex.axis=0.7,cex.lab=1,
  border=border,col=col3d,main="Exposure-lag-response")
mtext("GLM with a priori selection",cex=0.6)

par(mar=c(5,4,4,1))

plot(predslglm1,"overall",ylab="RR",xlab="Temperature (C)",xlim=c(-5,30),
  ylim=c(0.5,3.5),lwd=1.5,col=col,main="Overall cumulative exposure-response")
mtext("GLM with a priori selection",cex=0.6)

plot(predslglm1,var=29,xlab="Lag (days)",ylab="RR",ylim=c(0.9,1.4),
  lwd=1.5,col=col,main="Lag-response at 29C")
mtext("GLM with a priori selection",cex=0.6)

################################################################################
# GLM WITH AIC-BASED KNOT SELECTION

par(mar=c(1,1,3,1))

plot(pred3dglm2,xlab="Temperature (C)",zlab="RR",zlim=c(0.9,1.4),xlim=c(-5,30),
  ltheta=170,phi=35,lphi=30,shade=0.75,cex.axis=0.7,cex.lab=1,
  border=border,col=col3d,main="Exposure-lag-response")
mtext("GLM with AIC-based selection",cex=0.6)

par(mar=c(5,4,4,1))

plot(predslglm2,"overall",ylab="RR",xlab="Temperature (C)",xlim=c(-5,30),
  ylim=c(0.5,3.5),lwd=1.5,col=col,main="Overall cumulative exposure-response")
mtext("GLM with AIC-based selection",cex=0.6)

plot(predslglm2,var=29,xlab="Lag (days)",ylab="RR",ylim=c(0.9,1.4),
  lwd=1.5,col=col,main="Lag-response at 29C")
mtext("GLM with AIC-based selection",cex=0.6)

################################################################################
# GAM WITH ONE DIMENSION UNPENALIZED BY USING PARAMETRIC FUNCTIONS

par(mar=c(1,1,3,1))

plot(pred3dgam3,xlab="Temperature (C)",zlab="RR",zlim=c(0.9,1.4),xlim=c(-5,30),
  ltheta=170,phi=35,lphi=30,shade=0.75,cex.axis=0.7,cex.lab=1,
  border=border,col=col3d,main="Exposure-lag-response")
mtext("GAM with partial penalization",cex=0.6)

par(mar=c(5,4,4,1))

plot(predslgam3,"overall",ylab="RR",xlab="Temperature (C)",xlim=c(-5,30),
  ylim=c(0.5,3.5),lwd=1.5,col=col,main="Overall cumulative exposure-response")
mtext("GAM with partial penalization",cex=0.6)

plot(predslgam3,var=29,xlab="Lag (days)",ylab="RR",ylim=c(0.9,1.4),
  lwd=1.5,col=col,main="Lag-response at 29C")
mtext("GAM with partial penalization",cex=0.6)

dev.off()

#
