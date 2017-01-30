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
# SECOND EXAMPLE
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

pdf("fig4.pdf",height=6,width=9)
layout(matrix(1:6,ncol=3,byrow=TRUE))

################################################################################
# GAM

par(mar=c(1,1,3,1))

plot(predgam,xlab="WLM/year",ylab="Lag (years)",zlab="RR",zlim=c(0.95,1.3),
  theta=230,ltheta=230,lphi=45,cex.axis=0.7,cex.lab=0.8,nticks=9,
  border=border,col=col3d,main="Exposure-lag-response")
mtext("GAM with additional varying ridge penalty",cex=0.6)

par(mar=c(5,4,4,1))

plot(predgam,lag=15,xlab="WLM/year",ylab="RR",ylim=c(0.9,1.3),axes=T,lwd=1.5,
  col=col,main="Exposure-response at lag 15")
mtext("GAM with additional varying ridge penalty",cex=0.6)

plot(predgam,var=100,ylab="RR",xlab="Lag (years)",xlim=c(0,40),
  ylim=c(0.9,1.3),lwd=1.5,col=col,main="Lag-response at 100 WLM/Year")
mtext("GAM with additional varying ridge penalty",cex=0.6)

################################################################################
# GLM (AS IN GASPARRINI STAT MED 2014)

par(mar=c(1,1,3,1))

plot(predglm,xlab="WLM/year",ylab="Lag (years)",zlab="RR",zlim=c(0.95,1.3),
  theta=230,ltheta=230,lphi=45,cex.axis=0.7,cex.lab=0.8,nticks=9,
  border=border,col=col3d,main="Exposure-lag-response")
mtext("GLM with AIC-based selection",cex=0.6)

par(mar=c(5,4,4,1))

plot(predglm,lag=15,xlab="WLM/year",ylab="RR",ylim=c(0.9,1.3),axes=T,lwd=1.5,
  col=col,main="Exposure-response at lag 15")
mtext("GLM with AIC-based selection",cex=0.6)

plot(predglm,var=100,ylab="RR",xlab="Lag (years)",xlim=c(0,40),
  ylim=c(0.9,1.3),lwd=1.5,col=col,main="Lag-response at 100 WLM/Year")
mtext("GLM with AIC-based selection",cex=0.6)


dev.off()

#
