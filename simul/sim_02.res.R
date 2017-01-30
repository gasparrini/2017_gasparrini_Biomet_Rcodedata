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
# SIMULATION STUDY
# SUMMARIZE THE RESULTS IN TABLE AND GRAPHS
################################################################################

################################################################################
# TABLE 1

tab <- formatC(cbind(
  colMeans(sapply(timemodel,colMeans)),
  do.call(cbind,lapply(seq(combsim), function(j) 
    t(sapply(seq(bias), function(i) c(
      mean(edf[[i]][,j]),
      #mean(abs(bias[[i]][[j]]))/mean(abs(bias[[1]][[j]])),
      mean(cov[[i]][[j]][seq(0,10,0.25)!=cen[j]]),
      mean(rmse[[i]][[j]])/mean(rmse[[1]][[j]])
    )))
  ))
),format="f",digits=2)

colnames(tab) <- c('time',t(outer(paste0('s',1:3),
  c('edf','cov','rmse'),paste,sep='-')))

library(xtable)
xtable(tab)

# ANALYSIS OF CONVERGENCE
lapply(nonconv,colMeans)

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
# PLOT OF THE FUNCTIONS USED TO SIMULATE THE EFFECT SURFACE

pdf("figS1.pdf",height=4,width=12)

layout(matrix(1:2,ncol=2))
par(mar=c(5,4,1,1))

plot(0:100/10,flin(0:100/10)-flin(2),type="l",xlab="x",ylab=expression(f(x)),
  xlim=c(-0.5,10.5),ylim=c(-0.2,1),col=2,lty=2,lwd=2)
lines(0:100/10,fflex(0:100/10)-fflex(5),col=3,lty=4,lwd=2)
lines(0:100/10,fdnorm(0:100/10)-fdnorm(5),col=4,lty=5,lwd=2)
abline(h=0)
legend("topleft",c("flin","fflex","fdnorm"),lty=c(2,4,5),lwd=2,
  col=2:4,cex=0.8,bty="n")

# PLOT OF THE FUNCTIONS USED TO SIMULATE THE LAG STRUCTURE
plot(0:40,wconst(0:40),type="l",xlab="Lag",ylab=expression(w(l)),
  xlim=c(-0.5,40.5),ylim=c(-0.2,1),col=2,lty=2,lwd=2)
lines(0:40,wdecay(0:40),col=3,lty=4,lwd=2)
lines(0:40,wpeak1(0:40),col=4,lty=5,lwd=2)
lines(0:40,wpeak2(0:40),col=6,lty=6,lwd=2)
lines(0:40,wdnorm(0:40),col=7,lty=7,lwd=2)
legend("topright",c("wconst","wdecay","wpeak1","wpeak2","wdnorm"),
  lty=c(2,4:7),lwd=2,col=c(2:4,6,7),cex=0.8,bty="n")

dev.off()

################################################################################
# GRAPHICAL REPRESENTATION OF THE SIMULATED EFFECT SURFACES

# ARGUMENTS FOR 3D PLOTS
arg3D <- list(x=seq(0,10,0.25),y=0:40,ticktype="detailed",theta=230,
  ltheta=200,phi=30,lphi=30,xlab="x",ylab="Lag",zlab="log-RR",
  shade = 0.75,r=sqrt(3),d=5,cex.axis=0.7,cex.lab=0.8,border=border,
  col=col3d)

varval3d <- c(5,1,7)
lagval3d <- c(20,1,26)
varlines3d <- list(0:40/4,19:40/4,0:30/4)
laglines3d <- list(0:40,8:40,0:40)

zlimlow <- c(-0.007,-0.01,-0.003)
zlimhigh <- c(0.020,0.10,0.03)

pdf("fig1.pdf",height=3,width=9)

layout(matrix(1:3,ncol=3,byrow=TRUE))
par(mar=c(1,1,3,1))

for(j in seq(trueeff)) {
  d3 <- do.call(persp,modifyList(arg3D,list(z=trueeff[[j]],
    zlim=c(zlimlow[j],zlimhigh[j]))))
  title(names(combsim)[j],cex.main=1.5)
  lines (trans3d(x=varval3d[j],y=laglines3d[[j]],
    z=trueeff[[j]][as.character(varval3d[j]),paste0('lag',laglines3d[[j]])],
    pmat=d3),col=1,lwd=2)
  lines (trans3d(x=varlines3d[[j]],y=lagval3d[j],
    z=trueeff[[j]][as.character(varlines3d[[j]]),paste0('lag',lagval3d[[j]])],
    pmat=d3),col=1,lwd=2)
}

dev.off()

################################################################################
# GRAPH OF SIMULATION RESULTS

ylimllag <- c(-0.005,-0.02,-0.01)
ylimhlag <- c(0.015,0.04,0.03)
ylimlvar <- c(-0.006,-0.02,-0.01)
ylimhvar <- c(0.02,0.06,0.03)

pdf("fig2.pdf",height=12,width=8)
layout(matrix(1:21,7),heights=c(0.5,rep(1,6)))

for(j in seq(trueeff)) {
  
  par(mar=c(0,4,0,1))
  plot(1:10,type="n",frame.plot=F,axes=F,xlab="",ylab="")
  text(5,5,names(combsim)[j],cex=1.7,font=2)
  
  par(mar=c(5,4,0.5,1),mgp=c(2,1,0))
  
  for(i in c(1,2,6)) {
    plot(seq(0,10,0.25),trueeff[[j]][,paste0('lag',lagval3d[[j]])],type="n",
      xlab="x",ylab="log-RR",ylim=c(ylimlvar[j],ylimhvar[j]),frame.plot=FALSE)
    for(m in sample[[i]][[j]]) lines(seq(0,10,0.25),
      m[,paste0('lag',lagval3d[[j]])],lwd=1.5,col=grey(0.8))  
    abline(h=0)
    lines(seq(0,10,0.25),trueeff[[j]][,paste0('lag',lagval3d[[j]])],lwd=1.5)
    lines(seq(0,10,0.25),pred[[i]][[j]][,paste0('lag',lagval3d[[j]])],lty=5,
      col=col,lwd=1.5)
    mtext(modlab[[i]],cex=0.6)
  }
  
  for(i in c(1,2,6)) {
    plot(0:40,trueeff[[j]][as.character(varval3d[j]),],type="n",xlab="Lag",
      ylab="log-RR",ylim=c(ylimllag[j],ylimhlag[j]),frame.plot=FALSE)
    for(m in sample[[i]][[j]]) lines(0:40,m[as.character(varval3d[j]),],lwd=1.5,
      col=grey(0.8))
    abline(h=0)
    lines(0:40,trueeff[[j]][as.character(varval3d[j]),],lwd=1.5)
    lines(0:40,pred[[i]][[j]][as.character(varval3d[j]),],lty=5,col=col,lwd=1.5)
    mtext(modlab[[i]],cex=0.6)
  }
}

dev.off()

pdf("figS2.pdf",height=12,width=8)
layout(matrix(1:21,7),heights=c(0.5,rep(1,6)))

for(j in seq(trueeff)) {
  
  par(mar=c(0,4,0,1))
  plot(1:10,type="n",frame.plot=F,axes=F,xlab="",ylab="")
  text(5,5,names(combsim)[j],cex=1.7,font=2)
  
  par(mar=c(5,4,0.5,1),mgp=c(2,1,0))
  
  for(i in 3:5) {
    plot(seq(0,10,0.25),trueeff[[j]][,paste0('lag',lagval3d[[j]])],type="n",
      xlab="x",ylab="log-RR",ylim=c(ylimlvar[j],ylimhvar[j]),frame.plot=FALSE)
    for( m in sample[[i]][[j]]) lines(seq(0,10,0.25),
      m[,paste0('lag',lagval3d[[j]])],lwd=1.5,col=grey(0.8))  
    abline(h=0)
    lines(seq(0,10,0.25),trueeff[[j]][,paste0('lag',lagval3d[[j]])],lwd=1.5)
    lines(seq(0,10,0.25),pred[[i]][[j]][,paste0('lag',lagval3d[[j]])],lty=5,
      col=col,lwd=1.5)
    mtext(modlab[[i]],cex=0.6)
  }
  
  for(i in 3:5) {
    plot(0:40,trueeff[[j]][as.character(varval3d[j]),],type="n",xlab="Lag",
      ylab="log-RR",ylim=c(ylimllag[j],ylimhlag[j]),frame.plot=FALSE)
    for(m in sample[[i]][[j]]) lines(0:40,m[as.character(varval3d[j]),],lwd=1.5,
      col=grey(0.8))
    abline(h=0)
    lines(0:40,trueeff[[j]][as.character(varval3d[j]),],lwd=1.5)
    lines(0:40,pred[[i]][[j]][as.character(varval3d[j]),],lty=5,col=col,lwd=1.5)
    mtext(modlab[[i]],cex=0.6)
  }
}

dev.off()

pdf("figS3.pdf",height=8.5,width=8)
layout(matrix(1:15,5),heights=c(0.5,rep(1,4)))

for(j in seq(trueeff)) {
  
  par(mar=c(0,4,0,1))
  plot(1:10,type="n",frame.plot=F,axes=F,xlab="",ylab="")
  text(5,5,names(combsim)[j],cex=1.7,font=2)
  
  par(mar=c(5,4,0.5,1),mgp=c(2,1,0))
  
  for(i in 7:8) {
    plot(seq(0,10,0.25),trueeff[[j]][,paste0('lag',lagval3d[[j]])],type="n",
      xlab="x",ylab="log-RR",ylim=c(ylimlvar[j],ylimhvar[j]),frame.plot=FALSE)
    for( m in sample[[i]][[j]]) lines(seq(0,10,0.25),
      m[,paste0('lag',lagval3d[[j]])],lwd=1.5,col=grey(0.8))  
    abline(h=0)
    lines(seq(0,10,0.25),trueeff[[j]][,paste0('lag',lagval3d[[j]])],lwd=1.5)
    lines(seq(0,10,0.25),pred[[i]][[j]][,paste0('lag',lagval3d[[j]])],lty=5,
      col=col,lwd=1.5)
    mtext(modlab[[i]],cex=0.6)
  }
  
  for(i in 7:8) {
    plot(0:40,trueeff[[j]][as.character(varval3d[j]),],type="n",xlab="Lag",
      ylab="log-RR",ylim=c(ylimllag[j],ylimhlag[j]),frame.plot=FALSE)
    for(m in sample[[i]][[j]]) lines(0:40,m[as.character(varval3d[j]),],lwd=1.5,
      col=grey(0.8))
    abline(h=0)
    lines(0:40,trueeff[[j]][as.character(varval3d[j]),],lwd=1.5)
    lines(0:40,pred[[i]][[j]][as.character(varval3d[j]),],lty=5,col=col,lwd=1.5)
    mtext(modlab[[i]],cex=0.6)
  }
}

dev.off()

################################################################################
# GRAPH OF COVERAGE IN THE BI-DIMENSIONAL SURFACE

col <- colorRampPalette(c("black","white"))

for(i in c(seq(cov))) {
  
  for(j in c(seq(cov[[i]]))) {
    
    pdf(paste(names(cov[[i]])[j],names(cov)[i],".pdf",sep=""),width=5,height=4)
    
    par(mar=c(4,4,3,0.5))
    
    filled.contour(seq(0,10,0.25),0:40,cov[[i]][[j]],zlim=c(0,1),color.palette=col,
      xlab="Exposure",ylab="Lag",key.axes=axis(4,at=0:10/10),nlevels=10)
    title(names(cov[[i]])[j])
    mtext(modlab[[i]],cex=0.8)
    
    dev.off()
  }
}

#
