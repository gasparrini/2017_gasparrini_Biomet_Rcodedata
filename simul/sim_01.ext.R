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
# PERFORM THE SIMULATIONS (EXTERNAL METHOD)
################################################################################

################################################################################
# START THE ITERATIONS
# 2 LEVELS:
#   - TYPE OF SURFACE (COMBINATIONS OF FUNCTIONS USED FOR SIMULATE EFFECTS)
#   - RANDOM GENERATED RESPONSE
# THEN FIT DIFFERENT TYPES OF MODELS

time <- proc.time()
# LOOP ACROSS SIMULATED EXPOSURE-LAG-RESPONSES
#for(j in seq(nrow(combsim))) {
for(j in seq(combsim)) {
  
  # PRINT
  cat("\n\n ",names(combsim)[j],"\n")
  
################################################################################
# LOOP ACROSS RANDOMLY SIMULATED DATA
  
  for(i in seq(nsim)) {
    
    # PRINT
    cat(i,"")
    
    # SET THE SEED
    seed <- 13041975 + i
    set.seed(seed)
    
    # SIMULATE THE DATA
    cumeff <- apply(Q,1,fcumeff,0:40,combsim[j])
    suppressWarnings(y <- rpois(length(x),exp((log(base[j])+cumeff))))
    
################################################################################
# FIT GAM-REML
    
    p <- 1
    
    # RUN THE MODEL
    mtime <- proc.time()
    cb <- crossbasis(Q[,1],lag=c(0,40),argvar=list(fun='ps',df=9),
      arglag=list(fun='ps'))
    cbgamPen <- cbPen(cb)
    model <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method='REML')
    # STORE THE TIME
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # MEASURE SIGNAL-TO-NOISE
    #cor(na.omit(log(y+1)),log(predict(model,type="response")+1))
        
    # PREDICT
    cp <- crosspred(cb,model,from=0,to=10,by=0.25,bylag=1,cen=cen[j])
    
    # STORE THE RESULTS
    pred[[p]][[j]] <- pred[[p]][[j]] + cp$matfit
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$matfit - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$matfit-qn*cp$matse & 
        trueeff[[j]] <= cp$matfit+qn*cp$matse)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$matfit - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$matfit))
    edf[[p]][i,j] <- sum(model$edf[-1])
    nonconv[[p]][i,j] <- !model$converged
    
################################################################################
# FIT GLM-AIC
    
    p <- 2
    
    # LOOP ACROSS GLM WITH DIFFERENT NUMBER OF KNOTS
    modellist <- list()
    mtime <- proc.time()
    for(k in seq(argvarlist1)) {
      
      # DEFINE THE CROSS-BASIS
      # NB: USE FIRST COLUMN OF Q AS TIME SERIES DATA -> FASTER
      cb <- crossbasis(Q[,1],lag=c(0,40),argvar=argvarlist1[[k]],
        arglag=arglaglist1[[k]])
      
      # RUN THE MODEL, SAVING IT IN THE LIST WITH MINIMAL INFO (SAVE MEMORY)
      modellist[[k]] <- glm(y~cb,family=poisson,model=F)
    }
    
    # DETERMINE THE BEST GLM FOR AIC
    best <- which.min(sapply(modellist,AIC))
    cb <- crossbasis(Q[,1],lag=c(0,40),argvar=argvarlist1[[best]],
      arglag=arglaglist1[[best]])
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # PREDICT
    cp <- crosspred(cb,modellist[[best]],from=0,to=10,by=0.25,bylag=1,cen=cen[j])
    
    # STORE THE RESULTS
    pred[[p]][[j]] <- pred[[p]][[j]] + cp$matfit
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$matfit - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$matfit-qn*cp$matse & 
        trueeff[[j]] <= cp$matfit+qn*cp$matse)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$matfit - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$matfit))
    edf[[p]][i,j] <- length(cp$coef)
    nonconv[[p]][i,j] <- !modellist[[k]]$converged
    
################################################################################
# FIT GAM-AIC
    
    p <- 3
    
    # RUN THE MODEL
    mtime <- proc.time()
    cb <- crossbasis(Q[,1],lag=c(0,40),argvar=list(fun='ps',df=9),
      arglag=list(fun='ps'))
    cbgamPen <- cbPen(cb)
    model <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method="GCV.Cp")
    # STORE THE TIME
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # PREDICT
    cp <- crosspred(cb,model,from=0,to=10,by=0.25,bylag=1,cen=cen[j])
    
    # STORE THE RESULTS
    pred[[p]][[j]] <- pred[[p]][[j]] + cp$matfit
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$matfit - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$matfit-qn*cp$matse & 
        trueeff[[j]] <= cp$matfit+qn*cp$matse)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$matfit - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$matfit))
    edf[[p]][i,j] <- sum(model$edf[-1])
    nonconv[[p]][i,j] <- !model$converged

################################################################################
# FIT GAMcr
    
    p <- 4
    
    # RUN THE MODEL
    mtime <- proc.time()
    cb <- crossbasis(Q[,1],lag=c(0,40),argvar=list(fun='cr',df=9),
      arglag=list(fun='cr'))
    cbgamPen <- cbPen(cb)
    model <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method='REML')
    # STORE THE TIME
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # PREDICT
    cp <- crosspred(cb,model,from=0,to=10,by=0.25,bylag=1,cen=cen[j])
    
    # STORE THE RESULTS
    pred[[p]][[j]] <- pred[[p]][[j]] + cp$matfit
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$matfit - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$matfit-qn*cp$matse & 
        trueeff[[j]] <= cp$matfit+qn*cp$matse)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$matfit - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$matfit))
    edf[[p]][i,j] <- sum(model$edf[-1])
    nonconv[[p]][i,j] <- !model$converged
    
################################################################################
# FIT GAMps21
    
    p <- 5
    
    # RUN THE MODEL
    mtime <- proc.time()
    cb <- crossbasis(Q[,1],lag=c(0,40),argvar=list(fun='ps',df=9),
      arglag=list(fun='ps',diff=1))
    cbgamPen <- cbPen(cb)
    model <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method='REML')
    # STORE THE TIME
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # PREDICT
    cp <- crosspred(cb,model,from=0,to=10,by=0.25,bylag=1,cen=cen[j])
    
    # STORE THE RESULTS
    pred[[p]][[j]] <- pred[[p]][[j]] + cp$matfit
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$matfit - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$matfit-qn*cp$matse & 
        trueeff[[j]] <= cp$matfit+qn*cp$matse)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$matfit - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$matfit))
    edf[[p]][i,j] <- sum(model$edf[-1])
    nonconv[[p]][i,j] <- !model$converged

################################################################################
# FIT GAM-last
    
    p <- 6
    
    # RUN THE MODEL
    mtime <- proc.time()
    cb <- crossbasis(Q[,1],lag=c(0,40),argvar=list(fun='ps',df=9),
      arglag=list(fun='ps'))
    cbgamPen <- cbPen(cb,addSlag=rep(0:1,c(6,4)))
    model <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method='REML')
    # STORE THE TIME
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # PREDICT
    cp <- crosspred(cb,model,from=0,to=10,by=0.25,bylag=1,cen=cen[j])
    
    # STORE THE RESULTS
    pred[[p]][[j]] <- pred[[p]][[j]] + cp$matfit
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$matfit - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$matfit-qn*cp$matse & 
        trueeff[[j]] <= cp$matfit+qn*cp$matse)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$matfit - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$matfit))
    edf[[p]][i,j] <- sum(model$edf[-1])
    nonconv[[p]][i,j] <- !model$converged
    
################################################################################
# FIT GAM-quad
    
    p <- 7
    
    # RUN THE MODEL
    mtime <- proc.time()
    cb <- crossbasis(Q[,1],lag=c(0,40),argvar=list(fun='ps',df=9),
      arglag=list(fun='ps',fx=T))
    C <- do.call(onebasis,c(list(x=0:40),attr(cb,"arglag")))
    D <- diff(diag(41),diff=2)
    P <- diag((0:38)^2)
    Slag2 <- t(C)%*%t(D)%*%P%*%D%*%C
    cbgamPen <- cbPen(cb,addSlag=Slag2)
    model <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method='REML')
    # STORE THE TIME
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # PREDICT
    cp <- crosspred(cb,model,from=0,to=10,by=0.25,bylag=1,cen=cen[j])
    
    # STORE THE RESULTS
    pred[[p]][[j]] <- pred[[p]][[j]] + cp$matfit
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$matfit - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$matfit-qn*cp$matse & 
        trueeff[[j]] <= cp$matfit+qn*cp$matse)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$matfit - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$matfit))
    edf[[p]][i,j] <- sum(model$edf[-1])
    nonconv[[p]][i,j] <- !model$converged

################################################################################
# FIT GAM-exp
  
    p <- 8
    
    # RUN THE MODEL
    mtime <- proc.time()
    cb <- crossbasis(Q[,1],lag=c(0,40),argvar=list(fun='ps',df=9),
      arglag=list(fun='ps',fx=T))
    cbgamPen <- cbPen(cb,addSlag=exp(0:9))
    model <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method='REML')
    # STORE THE TIME
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # PREDICT
    cp <- crosspred(cb,model,from=0,to=10,by=0.25,bylag=1,cen=cen[j])
    
    # STORE THE RESULTS
    pred[[p]][[j]] <- pred[[p]][[j]] + cp$matfit
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$matfit - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$matfit-qn*cp$matse & 
        trueeff[[j]] <= cp$matfit+qn*cp$matse)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$matfit - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$matfit))
    edf[[p]][i,j] <- sum(model$edf[-1])
    nonconv[[p]][i,j] <- !model$converged
    
  }
  
################################################################################
# COMPUTE THE AVERAGE FOR THE STATISTICS OF THE WHOLE SURFACE
  
  for(i in seq(pred)) {
    pred[[i]][[j]] <- pred[[i]][[j]]/nsim
    bias[[i]][[j]] <- bias[[i]][[j]]/nsim
    cov[[i]][[j]] <- cov[[i]][[j]]/nsim
    rmse[[i]][[j]] <- sqrt(rmse[[i]][[j]]/nsim)
  }

}
(tottime <- proc.time()-time)

################################################################################
# SAVE

save.image("simext.RData")

#
