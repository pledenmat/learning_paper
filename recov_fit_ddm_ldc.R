rm(list=ls())
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
library(myPackage)
library(gamlss.dist) # Ex-gaussian for generating confRT distribution
library(DEoptim)
source("ldc_nn_functions.R")
library(Rcpp)
sourceCpp("ldc_train.cpp")

ldc_conf <- function(a,b,data){
  return(( (1+data$resp)/2)*(1 /(1 + exp( (1/sqrt(data$rt2+.00001))* (-a*data$evidence2 - b))))  + ((1-data$resp)/2)*(1 /(1 + exp((1/sqrt(data$rt2+.00001))* (a*data$evidence2 - b))) ))
}

# Generating parameters ---------------------------------------------------
v <- .1 #' Drift rate
s <- .1 #' Within-trial noise
a <- .07 #' Boundary distance from 0
ter <- .4 #' Non-decision time
z <- 0 #' Starting point bias
dt <- .001 #' Time step
ntrials <- 200
t2time <- rexGAUS(ntrials,mu=.2,sigma=.1,nu=.3) #' Post-decision accumulation time
t2time[t2time<0] <- .1
vratio <- 1

alpha <- 15
beta <- 1
# Generate dataset --------------------------------------------------------

Data <- data.frame(DDM_with_confidence_slow_fullconfRT(v,a,ter,z,ntrials,s,dt,t2time,vratio))
names(Data) <- c('rt','resp','cor','evidence2','rt2','cj')
Data$difflevel <- 'medium'
Data$RTconf <- Data$rt2 - Data$rt
Data$evidence2 <- Data$evidence2*Data$resp
Data$resp <- 1

Data$cj <- ldc_conf(alpha,beta,Data)


# Fit dataset -------------------------------------------------------------

w0 <- c(10,0,1)
beta_input <- .1

file_name <- paste0('recovery/testfit.Rdata')
params_lower <- c(0,0,0,1,0)
params_upper <- c(.2,2,0,1,.5)
if(file.exists(file_name)){
  load(file_name)
}else{
  optimal_params <- DEoptim(ldc.nn.fit, # function to optimize
                            lower = params_lower, upper = params_upper,obs = Data,
                            dt = dt, sigma = s, Nupdate_per_trial = 1000,w0=w0,beta_input=beta_input,
                            control=c(itermax=1000,steptol=100,reltol=.001,NP=40))
  results <- summary(optimal_params)
  #save individual results
  save(results, file=file_name)
}


# For debugging/testing the functions -------------------------------------

# a,ter,z,vratio,v
w0 <- c(0,0,1)
ddm_params <- c(a,v)
ntrials=10;dt=.001;sigma=.1;Nupdate_per_trial=1000;
confRTname="RTconf";diffname="difflevel";respname="resp";
totRTname='rt2';targetname='cj';accname='cor';beta_input=.1;error_type='mse'
obs <- temp_dat
binning=FALSE;nbin=6
params <- ldc.results$optim$bestmem[1]
