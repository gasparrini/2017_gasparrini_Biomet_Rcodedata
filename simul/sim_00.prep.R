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
# PREPARE THE DATA
################################################################################

# LOAD PACKAGES AND FUNCTIONS
library(dlnm) ; library(mgcv) ; library(splines) ; library(tsModel)

################################################################################
# DEFINE THE EXPOSURE

# REAL TEMPERATURE SERIES, STANDARDIZED IN 0-10
x <- chicagoNMMAPS$temp
x <- (x-min(x))/diff(range(x))*10

# MATRIX Q OF EXPOSURE HISTORIES
Q <- Lag(x,0:40)

################################################################################
# DEFINE THE BI-DIMENSIONAL ASSOCIATIONS USED FOR SIMULATING DATA

# BASIC FUNCTIONS TO SIMULATE UNIDIMENSIONAL SHAPES
flin <- function(x) 0.1*x
fflex <- function(x) {
  coef <- c(0.2118881,0.1406585,-0.0982663,0.0153671,-0.0006265)
  as.numeric(outer(x,0:4,'^')%*%coef)
}
fdnorm <- function(x) (dnorm(x,1.5,2)+1.5*dnorm(x,7.5,1))

wconst <- function(lag) lag-lag+0.20
wdecay <- function(lag) exp(-lag/2)
wpeak1 <- function(lag) 12*dnorm(lag,8,5)
wpeak2 <- function(lag) 15*dnorm(lag,8,10)
wdnorm <- function(lag) 5*(dnorm(lag,4,6)+dnorm(lag,25,4))

# FUNCTIONS TO SIMULATE THE BI-DIMENSIONAL EXPOSURE-LAG-RESPONSE
fplane <- function(x,lag) 0.1 * (flin(x)-flin(2)) * wconst(lag)
ftemp <- function(x,lag) 0.1 * (fflex(x)-fflex(5)) * 
    ifelse(is.na(x),NA,ifelse(x>=5,wdecay(lag),wpeak1(lag)))
fcomplex <- function(x,lag) 0.1 * (fdnorm(x)-fdnorm(5)) * 
    ifelse(is.na(x),NA,ifelse(x>=5,wdnorm(lag),wpeak2(lag)))

# COMBINATIONS OF FUNCTIONS USED TO SIMULATE DATA
combsim <- c("fplane","ftemp","fcomplex")
names(combsim) <- c("Plane","Temperature","Complex")

# CENTERING POINT
cen <- c(2,5,5)

# LIST WITH TRUE EFFECT SURFACES OF FOR EACH COMBINATIONS
trueeff <- lapply(combsim, function(fun) {
  temp <- outer(seq(0,10,0.25),0:40,fun)
  dimnames(temp) <- list(seq(0,10,0.25),paste("lag",0:40,sep=""))
  return(temp)
})
names(trueeff) <- names(combsim)

# FUNCTION TO COMPUTE THE CUMULATIVE EFFECT GIVEN AN EXPOSURE HISTORY
fcumeff <- function(hist,lag,fun) sum(do.call(fun,list(hist,lag)))

################################################################################
# DEFINE OPTIONS FOR GLM MODELS

# MAXIMUM NUMBER OF KNOTS AND DEGREE IN GLM
kmax <- 7
degree <- 2

# SET OF OPTIONS: INCREASING NUMBER OF KNOTS
argvarlist1 <- c(list(
  list(fun='lin'),list(fun='bs',degree=degree)),
  lapply(seq(kmax), function(k) {
    knots <- equalknots(0:10,nk=k,fun='bs',degree=degree)
    list(fun='bs',knots=knots,degree=degree)
  })
)
arglaglist1 <- c(list(
  list(fun='strata',int=T),list(fun='bs',degree=degree,int=T)),
  lapply(seq(kmax), function(k) {
    knots <- equalknots(0:40,nk=k,fun='bs',degree=degree)
    list(fun='bs',knots=knots,degree=degree,int=T)
  })
)
argvarlist1 <- rep(argvarlist1,kmax+2)
arglaglist1 <- rep(arglaglist1,each=kmax+2)


################################################################################
# DEFINE THE SETTING FOR THE SIMULATION

# NUMBER OF ITERATIONS (SET TO 1000 TO REPLICATE THE RESULTS)
nsim <- 10

# NUMBER OF SAMPLES USED AS EXAMPLES OF INDIVIDUAL ESTIMATES
nsample <- min(nsim,25)

# BASELINE
base <- c(15,150,15)

# NOMINAL VALUE
qn <- qnorm(0.975)

################################################################################
# CREATE THE OBJECT TO STORE RESULTS

# MODELS
models <- c('GAM','GLM-AIC','GAM-AIC','GAM-cr','GAM-ps21','GAM-last','GAM-quad',
  'GAM-exp')
modlab <-list(expression(plain(GAM)),expression(plain(GLM-AIC)),
  expression(plain(GAM-AIC)),expression(plain(GAM-CR)),
  expression(plain(GAM-PS)[plain("2,1")]),expression(plain(GAM-ADD)[plain(LAST)]),
  expression(plain(GAM-ALT)[plain(QUAD)]),expression(plain(GAM-ALT)[plain(EXP)]))

# LISTS FOR STORING THE AVERAGE PREDICTION, BIAS, COVERAGE RMSE AND SAMPLES
# FOR THE WHOLE SURFACE
pred <- lapply(models,function(i) lapply(seq(trueeff), function(j) 0))
sample <- lapply(models,function(i) vector('list',length(trueeff)))
names(pred) <- names(sample) <- models
for(i in seq(pred)) names(pred[[i]]) <- names(sample[[i]]) <- names(trueeff)
bias <- cov <- rmse <- pred

# LISTS OF MATRICES FOR STORING EDF, TIME, AND CONVERGENCE
nonconv <- lapply(models,function(x) 
  matrix(NA,nsim,length(trueeff),dimnames=list(seq(nsim),names(pred[[1]]))))
names(nonconv) <- models
edf <- timemodel <- nonconv

#
