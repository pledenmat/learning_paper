rm(list=ls())
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
library(myPackage) # Run devtools::install_github("pledenmat/myPackage") to install this custom package
library(Rcpp)
library(DEoptim)
library(ggplot2)
sourceCpp("ldc_train.cpp")
source("ldc_nn_functions.R")
set.seed(666)

ldc <- function(a,b,data){
  return( 1 /(1 + exp( (1/sqrt(data$rt2+.00001))* (-a*data$evidence2*data$resp - b))))
}

# Load data and setup simulation parameters ---------------------------------------------
Data <- read.csv("alternating_fb_mod_trim.csv")
# write.csv(ddm_par[,c("bound","drift","ter","sub","condition","difflevel","manip")],file="ddm_params.csv",row.names = F)
ddm_par <- read.csv("fit/alternating_fb/ddm_params.csv")
gen_par_alpha_learn <- read.csv("fit/alternating_fb/par_alpha_learn.csv")
gen_par_both_learn <- read.csv("fit/alternating_fb/par_both_learn.csv")
z <- 0
sigma <- .1
dt <- .001
difflevels <- sort(unique(Data$difflevel))
subs <- sort(unique(Data$sub))
Nsub <- length(subs)
alpha_minus <- 9
alpha_neutral <- 18
alpha_plus <- 36
beta_minus <- -1
beta_neutral <- 0
beta_plus <- 1
Nupdate_per_trial <- 1
binning <- F
#' First, we generate trials from DDM parameters
rm(simDat)
for (s in 1:Nsub) {
  print(s)
  for(level in difflevels){
    temp_par <- subset(ddm_par,sub==subs[s]&difflevel==level)[1,]
    temp_dat <- subset(Data,sub==subs[s]&difflevel==level)
    ntrial <- nrow(temp_dat)
    predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(
      v=temp_par$drift,a=temp_par$bound,ter=temp_par$ter,z=z,ntrials=ntrial,
      s=sigma,dt=dt,t2distribution=temp_dat$RTconf,postdriftmod=1))
    names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
    predictions$difflevel <- level
    predictions$sub <- subs[s]
    predictions$manip <- unique(temp_dat$manip)
    predictions$switch <- 1:nrow(predictions) %% 6
    predictions$trial <- -99 # Actual trial numbers are given at the subject level below
    if (unique(temp_dat$group)=='plus_first') {
      predictions$condition <- "minus"
      predictions[predictions$switch %in% c(0,2,4),"condition"] <- "plus"
    } else if (unique(temp_dat$group)=='minus_first') {
      predictions$condition <- "plus"
      predictions[predictions$switch %in% c(0,2,4),"condition"] <- "minus"
    } else {
      print(paste("No group detected for subject",subs[s]))
    }
    if(level==difflevels[1]&s==1){ simDat <- predictions
    }else{
      simDat <- rbind(simDat,predictions)
    }
  }
  #' We reorganize trial order to reproduce the actual experiment, i.e.:
  #' - Intermixed trial difficulty
  #' - Alternating feedback conditions
  temp <- simDat[simDat$sub==subs[s],]
  temp <- temp[order(temp$sub,temp$switch),]
  count_switch <- table(temp$switch)
  cum_switch <- cumsum(count_switch)
  shuffle <- c(sample(1:cum_switch[1]),sample((cum_switch[1]+1):cum_switch[2]),
               sample((cum_switch[2]+1):cum_switch[3]),sample((cum_switch[3]+1):cum_switch[4]),
               sample((cum_switch[4]+1):cum_switch[5]),sample((cum_switch[5]+1):cum_switch[6]))
  simDat[simDat$sub==subs[s],] <- temp[shuffle,]
  simDat[simDat$sub==subs[s],"trial"] <- 1:nrow(simDat[simDat$sub==subs[s],])
}

# Then we provide feedback from same parameters than the experiment
simDat$fb <- -99
simDat[simDat$manip=='alpha'&simDat$condition=='minus','fb'] <- ldc(alpha_minus,beta_neutral,simDat[simDat$manip=='alpha'&simDat$condition=='minus',])
simDat[simDat$manip=='beta'&simDat$condition=='minus','fb'] <- ldc(alpha_neutral,beta_minus,simDat[simDat$manip=='beta'&simDat$condition=='minus',])
simDat[simDat$manip=='alpha'&simDat$condition=='plus','fb'] <- ldc(alpha_plus,beta_neutral,simDat[simDat$manip=='alpha'&simDat$condition=='plus',])
simDat[simDat$manip=='beta'&simDat$condition=='plus','fb'] <- ldc(alpha_neutral,beta_plus,simDat[simDat$manip=='beta'&simDat$condition=='plus',])

simDat$evidence <- simDat$evidence2 # Variable naming for generating confidence

# Don't have all the fits from subject 98 so skip them for now
# simDat <- subset(simDat,sub != subs[98])
subs <- sort(unique(simDat$sub))
Nsub <- length(subs)


#' Finally, we compute confidence.
#' 2 simulated datasets are generated:
#' - alpha learning rate taken from alpha model and beta learning rate = 0
#' - alpha learning rate taken from alpha model and beta learning rate taken from both model
#' For both datasets, a0 and b0 retrieved from both model
#' TO DO ONCE NEW MODEL ARE FITTED:
#' - simulate 4 datasets, one for each combination of learning rate or not
simDat_alpha <- simDat
simDat_both <- simDat
for (s in 1:Nsub) {
  print(s)
  ddm_params <- subset(ddm_par,sub==subs[s])
  ddm_params <- c(mean(ddm_params$bound),
                  with(ddm_params,aggregate(drift,list(difflevel),mean))$x)
  temp_par_both <- subset(gen_par_both_learn,sub==subs[s])
  temp_par_alpha <- subset(gen_par_alpha_learn,sub==subs[s])
  temp_dat <- subset(simDat,sub==subs[s])
  eta_a <- mean(subset(gen_par_alpha_learn,sub==subs[s])$eta_a)
  results_alpha_learn <-
    ldc.nn.fit.w(params=c(mean(temp_par_alpha$a0),
                          mean(temp_par_alpha$b0),1,
                          mean(temp_par_alpha$eta_a),
                          0),
                 ddm_params = ddm_params,
                 obs=temp_dat,returnFit = F,eta_sep=T,estimate_evidence = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma,targetname = 'fb',fitname = 'cj')
  # results_both_learn <-
  #   ldc.nn.fit.w(params=c(mean(temp_par_both$a0),
  #                         mean(temp_par_both$b0),1,
  #                         mean(temp_par_alpha$eta_a),
  #                         mean(temp_par_both$eta_b)),
  #                ddm_params = ddm_params,
  #                obs=temp_dat,returnFit = F,eta_sep=T,estimate_evidence = F,
  #                Nupdate_per_trial=Nupdate_per_trial, binning = binning,
  #                dt = dt, sigma = sigma,targetname = 'fb',fitname = 'cj')
  simDat_alpha[simDat_alpha$sub==subs[s],'cj'] <- results_alpha_learn$pred
  # simDat_both[simDat_both$sub==subs[s],'cj'] <- results_both_learn$pred
}

# Save generated datasets
simDat_alpha$RTconf <- simDat_alpha$rt2 - simDat_alpha$rt 
write.csv(simDat_alpha,"simDat_alpha2.csv",row.names = F)
# write.csv(simDat_both,"simDat_both.csv",row.names = F)
# Fit models --------------------------------------------------------------
beta_input <- .1
Nupdate_per_trial <- 1
Nsim_err <- 1000
# bound, ter, z, vratio, drifts
dt <- .001; sigma <- .1
error_type1 <- "cross-entropy"
error_type2 <- "mse"
target <- "fb"
fitname <- 'cj'

par_beta_learn <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                             sub = subs, manip=NA, model = "beta",Npar = 3,family = "continuous_conf")
par_alpha_learn <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                              sub = subs, manip=NA, model = "alpha",Npar = 3,family = "continuous_conf")
par_both_learn <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                             sub = subs, manip=NA, model = "both",Npar = 4,family = "continuous_conf")
par_no_learn <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                           sub = subs,manip=NA, model = "no",Npar = 2,family = "continuous_conf")
par_bin <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                      sub = subs, manip=NA, model = "bin",Npar = 4,family = "binned_conf")
par_bin_mid <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                          sub = subs,manip=NA, model = "bin_mid",Npar = 4,family = "binned_conf")
par_ev_unknown <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                             sub = subs,manip=NA, model = "ev_unknown",Npar = 4,family = "ev_unknown")
par_mean_ev <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                             sub = subs,manip=NA, model = "mean_ev",Npar = 4,family = "ev_unknown")


simDat_alpha_bin <- simDat_alpha
simDat_alpha_bin$cj_continuous <- simDat_alpha_bin$cj
simDat_alpha_bin$cj_bin <- as.numeric(cut(simDat_alpha_bin$cj,breaks=seq(0,1,length.out=7),include.lowest = T))
simDat_alpha_bin$cj <- simDat_alpha_bin$cj_bin/6
for (s in 1:Nsub) {
  temp_dat <- subset(simDat_alpha_bin,sub==subs[s])
  print(paste("Running participant",s,"of",Nsub))
  
  ## Fit models with binned confidence
  ddm_file1 <- paste0('fit/alternating_fb/ddm/ddmfit_',subs[s],'.Rdata')
  if(file.exists(ddm_file1)){
    load(ddm_file1)
  }else{
    optimal_params <- DEoptim(ddm.fit, # function to optimize
                              lower = params_lower, upper = params_upper,obs = temp_dat,
                              dt = dt, sigma = sigma,ntrials=20,
                              control=c(itermax=1000,steptol=50,reltol=.001,NP=50))
    ddm.results <- summary(optimal_params)
    #save individual results
    save(ddm.results, file=ddm_file1)
  }
  ddm_params <- ddm.results$optim$bestmem[c(1,5:length(ddm.results$optim$bestmem))]
  
  
  
  ldc_file <- paste0('fit/alternating_fb/ldc_nn/recovery/binned_confidence/sim_alpha_learn/both_model/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
  }else{
    optimal_params <- DEoptim(ldc.nn.fit.w,ddm_params=ddm_params, 
                              obs = temp_dat,
                              lower = c(0,-100,1,0,0), 
                              upper = c(100,100,1,10000,10000),
                              Nupdate_per_trial=Nupdate_per_trial,
                              dt = dt, sigma = sigma,binning=binning,eta_sep=T,estimate_evidence=F,
                              Nsim_err=Nsim_err,error_type1=error_type1,fitname=fitname,
                              error_type2=error_type2,targetname=target,cost="separated",
                              control=c(itermax=1000,steptol=70,NP=40))
    ldc.results <- summary(optimal_params)
    #save individual results
    save(ldc.results, file=ldc_file)
  }
  par_bin[par_bin$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
  par_bin[par_bin$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
  par_bin[par_bin$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
  par_bin[par_bin$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
  par_bin[par_bin$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  
  
  temp_dat$cj <- temp_dat$cj - 1/12 # Try fixing bin values to center of their range
  ldc_file <- paste0('fit/alternating_fb/ldc_nn/recovery/binned_confidence/sim_alpha_learn/both_model_bin_mid/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
  }else{
    optimal_params <- DEoptim(ldc.nn.fit.w,ddm_params=ddm_params, 
                              obs = temp_dat,
                              lower = c(0,-100,1,0,0), 
                              upper = c(100,100,1,10000,10000),
                              Nupdate_per_trial=Nupdate_per_trial,
                              dt = dt, sigma = sigma,binning=binning,eta_sep=T,estimate_evidence=F,
                              Nsim_err=Nsim_err,error_type1=error_type1,fitname=fitname,
                              error_type2=error_type2,targetname=target,cost="separated",
                              control=c(itermax=1000,steptol=70,NP=40))
    ldc.results <- summary(optimal_params)
    #save individual results
    save(ldc.results, file=ldc_file)
  }
  par_bin_mid[par_bin_mid$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
  par_bin_mid[par_bin_mid$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
  par_bin_mid[par_bin_mid$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
  par_bin_mid[par_bin_mid$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
  par_bin_mid[par_bin_mid$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  
  
  ## Fit models with continuous confidence
  temp_dat <- subset(simDat_alpha,sub==subs[s])
  
  ldc_file <- paste0('fit/alternating_fb/ldc_nn/recovery/continuous_conf/sim_alpha_learn/beta_model/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
  }else{
    optimal_params <- DEoptim(ldc.nn.fit.w,ddm_params=ddm_params, 
                              obs = temp_dat,
                              lower = c(0,-100,1,0,0), 
                              upper = c(50,100,1,0,10000),
                              Nupdate_per_trial=Nupdate_per_trial,
                              dt = dt, sigma = sigma,binning=binning,eta_sep=T,estimate_evidence=F,
                              Nsim_err=Nsim_err,error_type1=error_type1,fitname=fitname,
                              error_type2=error_type2,targetname=target,cost="separated",
                              control=c(itermax=1000,steptol=70,NP=30))
    ldc.results <- summary(optimal_params)
    #save individual results
    save(ldc.results, file=ldc_file)
  }
  par_beta_learn[par_beta_learn$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
  par_beta_learn[par_beta_learn$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
  par_beta_learn[par_beta_learn$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
  par_beta_learn[par_beta_learn$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
  par_beta_learn[par_beta_learn$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  
  ldc_file <- paste0('fit/alternating_fb/ldc_nn/recovery/continuous_conf/sim_alpha_learn/alpha_model/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
  }else{
    optimal_params <- DEoptim(ldc.nn.fit.w,ddm_params=ddm_params, 
                              obs = temp_dat,
                              lower = c(0,-100,1,0,0), 
                              upper = c(50,100,1,10000,0),
                              Nupdate_per_trial=Nupdate_per_trial,
                              dt = dt, sigma = sigma,binning=binning,eta_sep=T,estimate_evidence=F,
                              Nsim_err=Nsim_err,error_type1=error_type1,fitname=fitname,
                              error_type2=error_type2,targetname=target,cost="separated",
                              control=c(itermax=1000,steptol=70,NP=30))
    ldc.results <- summary(optimal_params)
    #save individual results
    save(ldc.results, file=ldc_file)
  }
  par_alpha_learn[par_alpha_learn$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
  par_alpha_learn[par_alpha_learn$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
  par_alpha_learn[par_alpha_learn$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
  par_alpha_learn[par_alpha_learn$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
  par_alpha_learn[par_alpha_learn$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  
  ldc_file <- paste0('fit/alternating_fb/ldc_nn/recovery/continuous_conf/sim_alpha_learn/no_model/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
  }else{
    optimal_params <- DEoptim(ldc.nn.fit.w,ddm_params=ddm_params, 
                              obs = temp_dat,
                              lower = c(0,-100,1,0,0), 
                              upper = c(50,100,1,0,0),
                              Nupdate_per_trial=Nupdate_per_trial,
                              dt = dt, sigma = sigma,binning=binning,eta_sep=T,estimate_evidence=F,
                              Nsim_err=Nsim_err,error_type1=error_type1,fitname=fitname,
                              error_type2=error_type2,targetname=target,cost="separated",
                              control=c(itermax=1000,steptol=70,NP=30))
    ldc.results <- summary(optimal_params)
    #save individual results
    save(ldc.results, file=ldc_file)
  }
  par_no_learn[par_no_learn$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
  par_no_learn[par_no_learn$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
  par_no_learn[par_no_learn$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
  par_no_learn[par_no_learn$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
  par_no_learn[par_no_learn$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  
  ldc_file <- paste0('fit/alternating_fb/ldc_nn/recovery/continuous_conf/sim_alpha_learn/both_model/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
  }else{
    optimal_params <- DEoptim(ldc.nn.fit.w,ddm_params=ddm_params, 
                              obs = temp_dat,
                              lower = c(0,-100,1,0,0), 
                              upper = c(100,100,1,10000,10000),
                              Nupdate_per_trial=Nupdate_per_trial,
                              dt = dt, sigma = sigma,binning=binning,eta_sep=T,estimate_evidence=F,
                              Nsim_err=Nsim_err,error_type1=error_type1,fitname=fitname,
                              error_type2=error_type2,targetname=target,cost="separated",
                              control=c(itermax=1000,steptol=70,NP=40))
    ldc.results <- summary(optimal_params)
    #save individual results
    save(ldc.results, file=ldc_file)
  }
  par_both_learn[par_both_learn$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
  par_both_learn[par_both_learn$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
  par_both_learn[par_both_learn$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
  par_both_learn[par_both_learn$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
  par_both_learn[par_both_learn$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  
  ## Retrieve models with evidence unknown
  ldc_file_ev_unknown <- paste0('fit/alternating_fb/ldc_nn/recovery/evidence_unknown/sim_alpha_learn/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file_ev_unknown)) {
    load(ldc_file_ev_unknown)
    par_ev_unknown[par_ev_unknown$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
    par_ev_unknown[par_ev_unknown$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
    par_ev_unknown[par_ev_unknown$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
    par_ev_unknown[par_ev_unknown$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
    par_ev_unknown[par_ev_unknown$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
    par_ev_unknown[par_ev_unknown$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  }
  
  ldc_file_mean_ev <- paste0('fit/alternating_fb/ldc_nn/recovery/evidence_unknown/mean_ev/sim_alpha_learn/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file_mean_ev)) {
    load(ldc_file_mean_ev)
    par_mean_ev[par_mean_ev$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
    par_mean_ev[par_mean_ev$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
    par_mean_ev[par_mean_ev$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
    par_mean_ev[par_mean_ev$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
    par_mean_ev[par_mean_ev$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
    par_mean_ev[par_mean_ev$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  }
  
  par_bin[par_bin$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  par_bin_mid[par_bin_mid$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  par_both_learn[par_both_learn$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  par_beta_learn[par_beta_learn$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  par_alpha_learn[par_alpha_learn$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  par_no_learn[par_no_learn$sub==subs[s],"manip"] <- unique(temp_dat$manip)
}
par <- rbind(par_alpha_learn,par_beta_learn,par_both_learn,par_no_learn,
             par_bin,par_bin_mid,par_ev_unknown,par_mean_ev)
par$Ndata_point <- round(nrow(Data)/Nsub)

# Plot parameters ------------------------------------------

# Plot parameters and cost distributions for all models
ggplot(par, aes(x = a0, colour = model, label = model)) +
  geom_density() +
  theme_bw()
ggplot(par, aes(x = b0, colour = model, label = model)) +
  geom_density() +
  theme_bw()
ggplot(subset(par,!(model %in% c("beta","no"))), aes(x = eta_a, colour = model, label = model)) +
  geom_density() +
  theme_bw()
ggplot(subset(par,!(model %in% c("alpha","no"))), aes(x = eta_b, colour = model, label = model)) +
  geom_density() +
  theme_bw()
ggplot(subset(par,model %in% c("beta","ev_unknown")), aes(x = eta_b, colour = model, label = model)) +
  geom_density() +
  theme_bw()

# Plot parameters and cost distributions for binned confidence models only
ggplot(subset(par,family=="binned_conf"), aes(x = a0, colour = model, label = model)) +
  geom_density() +
  theme_bw()
ggplot(subset(par,family=="binned_conf"), aes(x = b0, colour = model, label = model)) +
  geom_density() +
  theme_bw()
ggplot(subset(par,family=="binned_conf"), aes(x = eta_a, colour = model, label = model)) +
  geom_density() +
  theme_bw()
ggplot(subset(par,family=="binned_conf"), aes(x = eta_b, colour = model, label = model)) +
  geom_density() +
  theme_bw()


# Formal model comparison -------------------------------------------------
bic_custom <- function(Residuals,k,n){
  return(log(n)*k+n*log(Residuals/n))
}

par$bic <- bic_custom(par$cost_ldc,par$Npar,par$Ndata_point)
mean_bic <- with(par,aggregate(bic,by=list(model=model,manip=manip),mean))
mean_bic$delta <- -99
mean_bic[mean_bic$manip=="alpha","delta"] <- 
  mean_bic[mean_bic$manip=="alpha",]$x - 
  min(mean_bic[mean_bic$manip=="alpha",]$x)
mean_bic[mean_bic$manip=="beta","delta"] <- 
  mean_bic[mean_bic$manip=="beta",]$x -
  min(mean_bic[mean_bic$manip=="beta",]$x)

bic <- with(mean_bic,aggregate(delta,list(model),mean))
names(bic) <- c('model','delta')
bic_order <- c('no','alpha','beta','both','bin','bin_mid','ev_unknown')
bic$model <- factor(bic$model, levels = bic_order)
bic <- bic[order(bic$model),]
barplot(bic[1:4,]$delta,names.arg = bic_order[1:4],xlab="Model (learning rate)", ylab = "BIC")
barplot(bic[c(4,5,7),]$delta,names.arg = c("continuous","binned","ev unknown"),
        xlab="Model (2 learning rates)", ylab = "BIC")
# Scatterplot generative/recovered parameters from selected models -----------------------------
gen_par_alpha_learn <- aggregate(cbind(a0,b0,eta_a,eta_b) ~ sub, data = gen_par_alpha_learn,mean)

# Best overall model
plot(par_alpha_learn$a0 ~ gen_par_alpha_learn$a0, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("a0, r =",round(cor(par_alpha_learn$a0, gen_par_alpha_learn$a0),3)))
abline(coef=c(0,1))

plot(par_alpha_learn$b0 ~ gen_par_alpha_learn$b0, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("b0, r =",round(cor(par_alpha_learn$b0, gen_par_alpha_learn$b0),3)))
abline(coef=c(0,1))

plot(par_alpha_learn$eta_a ~ gen_par_alpha_learn$eta_a, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("eta_a, r =",round(cor(par_alpha_learn$eta_a, gen_par_alpha_learn$eta_a),3)))
abline(coef=c(0,1))

### Comparison between both binned confidence models
plot(par_bin_mid$a0 ~ gen_par_alpha_learn$a0, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("a0, r =",round(cor(par_bin_mid$a0, gen_par_alpha_learn$a0),3)))
abline(coef=c(0,1))

plot(par_bin_mid$b0 ~ gen_par_alpha_learn$b0, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("b0, r =",round(cor(par_bin_mid$b0, gen_par_alpha_learn$b0),3)))
abline(coef=c(0,1))

plot(par_bin_mid$eta_a ~ gen_par_alpha_learn$eta_a, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("eta_a, r =",round(cor(par_bin_mid$eta_a, gen_par_alpha_learn$eta_a),3)))
abline(coef=c(0,1))

# Not much of a difference with the previous bin values used 
plot(par_bin$a0 ~ gen_par_alpha_learn$a0, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("a0, r =",round(cor(par_bin$a0, gen_par_alpha_learn$a0),3)))
abline(coef=c(0,1))

plot(par_bin$b0 ~ gen_par_alpha_learn$b0, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("b0, r =",round(cor(par_bin$b0, gen_par_alpha_learn$b0),3)))
abline(coef=c(0,1))

plot(par_bin$eta_a ~ gen_par_alpha_learn$eta_a, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("eta_a, r =",round(cor(par_bin$eta_a, gen_par_alpha_learn$eta_a),3)))
abline(coef=c(0,1))

# Evidence unknown
plot(par_ev_unknown$a0 ~ gen_par_alpha_learn$a0, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("a0, r =",round(cor(par_ev_unknown$a0, gen_par_alpha_learn$a0),3)))
abline(coef=c(0,1))

plot(par_ev_unknown$b0 ~ gen_par_alpha_learn$b0, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("b0, r =",round(cor(par_ev_unknown$b0, gen_par_alpha_learn$b0),3)))
abline(coef=c(0,1))

plot(par_ev_unknown$eta_a ~ gen_par_alpha_learn$eta_a, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("eta_a, r =",round(cor(par_ev_unknown$eta_a, gen_par_alpha_learn$eta_a),3)))
abline(coef=c(0,1))

plot(par_ev_unknown$eta_b ~ gen_par_alpha_learn$eta_b, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("eta_b, r =",round(cor(par_ev_unknown$eta_b, gen_par_alpha_learn$eta_b),3)))
abline(coef=c(0,1))

# Evidence unknown
plot(par_mean_ev$a0 ~ gen_par_alpha_learn$a0, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("a0, r =",round(cor(par_mean_ev$a0, gen_par_alpha_learn$a0),3)))
abline(coef=c(0,1))

plot(par_mean_ev$b0 ~ gen_par_alpha_learn$b0, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("b0, r =",round(cor(par_mean_ev$b0, gen_par_alpha_learn$b0),3)))
abline(coef=c(0,1))

plot(par_mean_ev$eta_a ~ gen_par_alpha_learn$eta_a, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("eta_a, r =",round(cor(par_mean_ev$eta_a, gen_par_alpha_learn$eta_a),3)))
abline(coef=c(0,1))

hist(par_mean_ev$eta_b,bty='n',main='',xlab='eta_b')
# Investigate the beta learning rate clusters -----------------------------
test <- subset(par_ev_unknown,eta_b < 100)
test_gen <- subset(gen_par_alpha_learn, sub %in% test$sub)
plot(test$a0 ~ test_gen$a0, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("a0, r =",round(cor(test$a0, test_gen$a0),3)))
abline(coef=c(0,1))

plot(test$b0 ~ test_gen$b0, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("b0, r =",round(cor(test$b0, test_gen$b0),3)))
abline(coef=c(0,1))

plot(test$eta_a ~ test_gen$eta_a, bty='n', xlab = "Generative", ylab = "Recovered",
     main = paste("eta_a, r =",round(cor(test$eta_a, test_gen$eta_a),3)))
abline(coef=c(0,1))

ddm_par$group <- "bad"
ddm_par[ddm_par$sub %in% test$sub,"group"] <- "good"

with(ddm_par,aggregate(drift,list(difflevel,group),mean))
with(ddm_par,aggregate(drift,list(difflevel,group),sd))

with(ddm_par,aggregate(bound,list(group),mean))
with(ddm_par,aggregate(bound,list(group),sd))

with(ddm_par,aggregate(ter,list(group),mean))
with(ddm_par,aggregate(ter,list(group),sd))

Data$group_model <- "bad"
Data[Data$sub %in% test$sub,"group_model"] <- "good"
m <- lmer(data=Data,rt ~ group + (1|sub))
anova(m)
m <- lmer(data=Data,cj ~ group + (1|sub))
anova(m)
m <- glmer(data=Data,cor ~ group + (1|sub),family = 'binomial')
library(car)
Anova(m)
