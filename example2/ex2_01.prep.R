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
# PREPARE THE DATA
################################################################################

# LOAD THE PACKAGES AND FUNCTIONS
library(survival) ; library(dlnm) ; library(splines) ; library(mgcv)

# LOAD THE DATA
uminers <- read.csv("uminers.csv",row.names=1)
nrow(uminers)

# ORDER THE DATASET BY RECORD ID
uminers <- uminers[order(uminers$record),]
rownames(uminers) <- uminers$record

################################################################################
# EXPAND THE DATASET BY 1-YEAR TIME PERIODS

# MAIN DATASET WITH OUTCOME, AGE AND CALENDAR YEAR
# - EXPAND THE DATA FOR EACH SUBJECT, CREATING SERIES OF AGE AND YEAR
# - CREATE THE OUTCOME SERIES, WHICH IS 1 ONLY IN THE LAST YEAR FOR NON-CENSORED
main <- do.call('rbind',lapply(seq(nrow(uminers)),function(i) {
  age <- seq(round(uminers[i,'agest']),round(uminers[i,'ageexit']))
  year <- 1000 + uminers[i,'byr'] + age
  id <- rep(uminers[i,'record'],length(age))
  out <- rep(c(0,uminers[i,'ind']),c(length(age)-1,1))
  return(data.frame(id,out,age,year))
}))

# FUNCTION TO RE-CREATE THE FULL EXPOSURE EXPERIENCE FROM 5-YEARS PERIODS, GIVEN:
#   - AGE AT START OF EXPOSURE (>5)
#   - AGE AT END OF EXPOSURE (<85)
#   - EXPOSURE LEVELS AT AGES 0-4, 5-9, ..., 85-89 (18 VALUES)
# STEPS:
#   - CREATE THE CUT-OFF POINTS
#   - CONVERT AGE AT START AND END IN YEARS WITH APPROPRIATE ROUNDING
#   - SUBSTITUTE THE APPROPRIATE CUT-OFF POINTS WITH AGE AT START AND END
#   - COMPUTE THE AVERAGE EXPOSURE FOR EACH YEAR (ADDING 10 MORE YEARS TO 99)
# NB: THE FUNCTION IS NOT VECTORIZED: ACCEPTS TWO SCALARS AND A VECTOR AS ARGS

fexpfull <- function(start,end,exp) {
  age <- 0:18*5
  start <- max(round(start-0.5),6)
  end <- min(round(end+0.5),84)
  age[(age-start)<0&(start-age)<5] <- start
  age[(age-end)>0&(age-end)<5] <- end
  expfull <- c(rep(as.numeric(exp)/diff(age),diff(age)),rep(0,10))
  names(expfull) <- 0:99
  return(expfull)
}

# LAGGED EXPOSURE MATRICES FOR RADON AND SMOKING
# - RECREATE THE EXPOSURE PROFILE FROM 5-Y PERIODS ACCOUNTING FOR START/END AGE
# - CREATE THE LAGGED EXPOSURES MATCHING AGE WITH PROFILE
# - NB: MATCHING START FROM AGE+1, AS FIRST OBS IN PROFILE IS AGE=0
Qx <- do.call('rbind',lapply(seq(nrow(uminers)),function(i) {
  profile <- fexpfull(uminers[i,'rdnstar'],uminers[i,'rendage'],
    as.numeric(uminers[i,17:34]))
  age <- main[main$id==i,"age"]
  hist <- exphist(profile,age,lag=c(2,40))
}))
Qz <- do.call('rbind',lapply(seq(nrow(uminers)),function(i) {
  profile <- fexpfull(uminers[i,'smkstar'],uminers[i,'sendage'],
    as.numeric(uminers[i,35:52]))
  age <- main[main$id==i,"age"]
  hist <- exphist(profile,age,lag=c(2,40))
}))

################################################################################
# GENERATE COVARIATES

cbz <- crossbasis(Qz,lag=c(2,40),argvar=list(fun="ns",knots=2.5),
  arglag=list(fun="strata",breaks=20))

#
