rm(list=ls())
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
library(myPackage) # Run devtools::install_github("pledenmat/myPackage") to install this custom package
library(Rcpp)
library(DEoptim)
library(ggplot2)
library(reshape)
sourceCpp("ldc_train.cpp")
source("ldc_nn_functions.R")
set.seed(666)

ldc <- function(a,b,data){
  return( 1 /(1 + exp( (1/sqrt(data$rt2+.00001))* (-a*data$evidence2*data$resp - b))))
}

# Retrieve estimated parameters from each model ---------------------------------------------
Data <- read.csv("alternating_fb_mod_trim.csv")

conditions <- sort(unique(Data$condition))
difflevels <- sort(unique(Data$difflevel)) # Important to sort to match with drift order
subs <- sort(unique(Data$sub))

totlen <- length(conditions)*length(difflevels)*length(subs)

# DDM parameters
ddm_par <- data.frame(bound = NA, drift = NA, ter = NA,cost_ddm = NA, 
                      sub = rep(subs, each = length(conditions)*length(difflevels)),
                      condition = rep(conditions,each = length(difflevels),
                                      length.out=totlen),
                      difflevel = rep(difflevels, length.out = totlen))

par_no_learn <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                           sub = rep(subs, each = length(conditions)*length(difflevels)),
                           condition = rep(conditions,each = length(difflevels),
                                           length.out=totlen),
                           difflevel = rep(difflevels, length.out = totlen),
                           manip=NA)

par_beta_learn <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                             sub = rep(subs, each = length(conditions)*length(difflevels)),
                             condition = rep(conditions,each = length(difflevels),
                                             length.out=totlen),
                             difflevel = rep(difflevels, length.out = totlen),
                             manip=NA)
par_alpha_learn <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                              sub = rep(subs, each = length(conditions)*length(difflevels)),
                              condition = rep(conditions,each = length(difflevels),
                                              length.out=totlen),
                              difflevel = rep(difflevels, length.out = totlen),
                              manip=NA)
par_both_learn <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                             sub = rep(subs, each = length(conditions)*length(difflevels)),
                             condition = rep(conditions,each = length(difflevels),
                                             length.out=totlen),
                             difflevel = rep(difflevels, length.out = totlen),
                             manip=NA)


go_to("fit")
go_to("alternating_fb")

for (s in 1:length(subs)) {
  print(paste("Retrieving participant",s))
  temp_dat <- subset(Data,sub==subs[s])
  
  ddm_file <- paste0('ddm/ddmfit_',subs[s],'.Rdata')
  if(file.exists(ddm_file)){
    load(ddm_file)
  }else{
    next
  }
  ddm_params <- ddm.results$optim$bestmem[c(1,5:length(ddm.results$optim$bestmem))]
  
  ddm_par[ddm_par$sub==subs[s],"bound"] <- ddm.results$optim$bestmem[1]
  ddm_par[ddm_par$sub==subs[s],"ter"] <- ddm.results$optim$bestmem[2]
  ddm_par[ddm_par$sub==subs[s],"z"] <- ddm.results$optim$bestmem[3]
  ddm_par[ddm_par$sub==subs[s],"vratio"] <- ddm.results$optim$bestmem[4]
  ddm_par[ddm_par$sub==subs[s]&ddm_par$difflevel==difflevels[1],"drift"] <- ddm.results$optim$bestmem[5]
  ddm_par[ddm_par$sub==subs[s]&ddm_par$difflevel==difflevels[2],"drift"] <- ddm.results$optim$bestmem[6]
  ddm_par[ddm_par$sub==subs[s]&ddm_par$difflevel==difflevels[3],"drift"] <- ddm.results$optim$bestmem[7]
  ddm_par[ddm_par$sub==subs[s],"cost_ddm"] <- ddm.results$optim$bestval
  
  ldc_file <- paste0('ldc_nn/trim/no_learn/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
    par_no_learn[par_no_learn$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
    par_no_learn[par_no_learn$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
    par_no_learn[par_no_learn$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
    par_no_learn[par_no_learn$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[4]
    par_no_learn[par_no_learn$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  }
  
  
  ldc_file_beta_learn <- paste0('ldc_nn/trim/beta_learn/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file_beta_learn)) {
    load(ldc_file_beta_learn)
    par_beta_learn[par_beta_learn$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
    par_beta_learn[par_beta_learn$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
    par_beta_learn[par_beta_learn$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
    par_beta_learn[par_beta_learn$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
    par_beta_learn[par_beta_learn$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  }
  
  ldc_file_alpha_learn <- paste0('ldc_nn/trim/alpha_learn/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file_alpha_learn)) {
    load(ldc_file_alpha_learn)
    par_alpha_learn[par_alpha_learn$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
    par_alpha_learn[par_alpha_learn$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
    par_alpha_learn[par_alpha_learn$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
    par_alpha_learn[par_alpha_learn$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
    par_alpha_learn[par_alpha_learn$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  }
  
  ldc_file_both_learn <- paste0('ldc_nn/trim/both_learn/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file_both_learn)) {
    load(ldc_file_both_learn)
    par_both_learn[par_both_learn$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
    par_both_learn[par_both_learn$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
    par_both_learn[par_both_learn$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
    par_both_learn[par_both_learn$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
    par_both_learn[par_both_learn$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  }
  
  par_no_learn[par_no_learn$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  par_both_learn[par_both_learn$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  par_beta_learn[par_beta_learn$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  par_alpha_learn[par_alpha_learn$sub==subs[s],"manip"] <- unique(temp_dat$manip)
}

# Merge parameters from every model
par_no_learn$model <- "no"
par_alpha_learn$model <- "alpha"
par_beta_learn$model <- "beta"
par_both_learn$model <- "both"

par <- rbind(par_alpha_learn,par_beta_learn,par_both_learn,par_no_learn)
par <- merge(par,ddm_par,by = c("sub","condition","difflevel"))


# Generate new datasets from these models & parameters -----------------------------

z <- 0
sigma <- .1
dt <- .001
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
    if (unique(temp_dat$order)=='plus_first') {
      predictions$condition <- "minus"
      predictions[predictions$switch %in% c(0,2,4),"condition"] <- "plus"
    } else if (unique(temp_dat$order)=='minus_first') {
      predictions$condition <- "plus"
      predictions[predictions$switch %in% c(0,2,4),"condition"] <- "minus"
    } else {
      print(paste("No order column detected for subject",subs[s]))
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

#' Finally, we compute confidence.
#' 4 simulated datasets are generated, one for each model.
#' Models differ in whether alpha/beta learning rate is fixed to 0 or not
simDat_alpha <- simDat
simDat_both <- simDat
simDat_beta <- simDat
simDat_no <- simDat
for (s in 1:Nsub) {
  print(s)
  ddm_params <- subset(ddm_par,sub==subs[s])
  ddm_params <- c(mean(ddm_params$bound),
                  with(ddm_params,aggregate(drift,list(difflevel),mean))$x)
  temp_par_both <- subset(par_both_learn,sub==subs[s])
  temp_par_alpha <- subset(par_alpha_learn,sub==subs[s])
  temp_par_beta <- subset(par_beta_learn,sub==subs[s])
  temp_par_no <- subset(par_no_learn,sub==subs[s])
  temp_dat <- subset(simDat,sub==subs[s])
  results_alpha_learn <-
    ldc.nn.fit.w(params=c(mean(temp_par_alpha$a0),
                          mean(temp_par_alpha$b0),1,
                          mean(temp_par_alpha$eta_a),
                          mean(temp_par_alpha$eta_b)),
                 ddm_params = ddm_params,
                 obs=temp_dat,returnFit = F,eta_sep=T,estimate_evidence = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma,targetname = 'fb',fitname = 'cj')
  results_both_learn <-
    ldc.nn.fit.w(params=c(mean(temp_par_both$a0),
                          mean(temp_par_both$b0),1,
                          mean(temp_par_both$eta_a),
                          mean(temp_par_both$eta_b)),
                 ddm_params = ddm_params,
                 obs=temp_dat,returnFit = F,eta_sep=T,estimate_evidence = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma,targetname = 'fb',fitname = 'cj')
  results_beta_learn <-
    ldc.nn.fit.w(params=c(mean(temp_par_beta$a0),
                          mean(temp_par_beta$b0),1,
                          mean(temp_par_beta$eta_a),
                          mean(temp_par_beta$eta_b)),
                 ddm_params = ddm_params,
                 obs=temp_dat,returnFit = F,eta_sep=T,estimate_evidence = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma,targetname = 'fb',fitname = 'cj')
  results_no_learn <-
    ldc.nn.fit.w(params=c(mean(temp_par_no$a0),
                          mean(temp_par_no$b0),1,
                          mean(temp_par_no$eta_a),
                          mean(temp_par_no$eta_b)),
                 ddm_params = ddm_params,
                 obs=temp_dat,returnFit = F,eta_sep=T,estimate_evidence = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma,targetname = 'fb',fitname = 'cj')
  
  simDat_alpha[simDat_alpha$sub==subs[s],'cj'] <- results_alpha_learn$pred
  simDat_both[simDat_both$sub==subs[s],'cj'] <- results_both_learn$pred
  simDat_no[simDat_no$sub==subs[s],'cj'] <- results_no_learn$pred
  simDat_beta[simDat_beta$sub==subs[s],'cj'] <- results_beta_learn$pred
}

# Save generated datasets
simDat_alpha$RTconf <- simDat_alpha$rt2 - simDat_alpha$rt 
simDat_no$RTconf <- simDat_no$rt2 - simDat_no$rt 
simDat_beta$RTconf <- simDat_beta$rt2 - simDat_beta$rt 
simDat_both$RTconf <- simDat_both$rt2 - simDat_both$rt 
# write.csv(simDat_alpha,"simDat_alpha.csv",row.names = F)
# write.csv(simDat_both,"simDat_both.csv",row.names = F)
# write.csv(simDat_beta,"simDat_beta.csv",row.names = F)
# write.csv(simDat_no,"simDat_no.csv",row.names = F)
# Retrieve model fits with evidence unknown --------------------------------------------------------------

models <- c("no","alpha","beta","both")

fit_par <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                      sub = rep(subs,each=length(models)^2), manip=NA, gen_model = rep(models,length(models)),
                      fit_model = rep(models,each = length(models)),Npar = NA,
                      Ndata_point = round(nrow(Data)/Nsub))
fit_par$Npar <- 3
fit_par[fit_par$fit_model=="no",'Npar'] <- 2
fit_par[fit_par$fit_model=="both",'Npar'] <- 4
for (s in 1:Nsub) {
  temp_dat <- subset(simDat,sub==subs[s])
  print(paste("Running participant",s,"of",Nsub))
  
  ## Fit models with binned confidence
  ddm_file1 <- paste0('ddm/ddmfit_',subs[s],'.Rdata')
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
  
  for (gen_model in models) {
    for (fit_model in models) {
      ldc_file <- paste0('ldc_nn/recovery/evidence_unknown/sim_',gen_model,'_learn/',fit_model,'_model/ldcfit_',subs[s],'.Rdata')
      if (file.exists(ldc_file)) {
        load(ldc_file)
      } else {
        # Some fits were not properly saved so retrieved from cluster job output logs
        ldc_save <- paste0('ldc_nn/recovery/evidence_unknown/sim_',gen_model,'_learn/',fit_model,'_model/slurm-55402977_',s,'.out')
        if (file.exists(ldc_save)) {
          test <- readLines(ldc_save)
          check_finish <- "***** summary of DEoptim object ***** "
          if (test[length(test)-5] == check_finish) {
            recov_par <- unlist(strsplit(test[length(test)-4],split=" "))
            a0 <- as.numeric(recov_par[7])
            b0 <- as.numeric(recov_par[8])
            eta_a <- as.numeric(recov_par[10])
            eta_b <- as.numeric(recov_par[11])
            cost <- as.numeric(unlist(strsplit(test[length(test)-3],split=" "))[8])
            ldc.results <- list(optim=list(bestmem=c(a0,b0,1,eta_a,eta_b),bestval=cost))
            save(ldc.results,file=ldc_file)
          }
        }
      }
      fit_par[fit_par$sub==subs[s] & fit_par$gen_model==gen_model & fit_par$fit_model==fit_model,"a0"] <- ldc.results$optim$bestmem[1]
      fit_par[fit_par$sub==subs[s] & fit_par$gen_model==gen_model & fit_par$fit_model==fit_model,"b0"] <- ldc.results$optim$bestmem[2]
      fit_par[fit_par$sub==subs[s] & fit_par$gen_model==gen_model & fit_par$fit_model==fit_model,"eta_a"] <- ldc.results$optim$bestmem[4]
      fit_par[fit_par$sub==subs[s] & fit_par$gen_model==gen_model & fit_par$fit_model==fit_model,"eta_b"] <- ldc.results$optim$bestmem[5]
      fit_par[fit_par$sub==subs[s] & fit_par$gen_model==gen_model & fit_par$fit_model==fit_model,"cost_ldc"] <- ldc.results$optim$bestval
      
    }
  }
  fit_par[fit_par$sub==subs[s],"manip"] <- unique(temp_dat$manip)
}

# Fit each dataset with all models, assuming trial evidence is known --------
beta_input <- .1
Nupdate_per_trial <- 1
Nsim_err <- 1000
# bound, ter, z, vratio, drifts
dt <- .001; sigma <- .1
error_type1 <- "cross-entropy"
error_type2 <- "mse"
target <- "fb"
fitname <- 'cj'
models <- c("no","alpha","beta","both")
fit_par_ev_known <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                               sub = rep(subs,each=length(models)^2), manip=NA, fit_model = models,
                               Npar = c(2,3,3,4), gen_model=rep(models,each=length(models)),
                               Ndata_point = round(nrow(Data)/Nsub))


simDat_alpha$gen_model <- "alpha" 
simDat_no$gen_model <- "no"
simDat_beta$gen_model <- "beta"
simDat_both$gen_model <- "both"

simDat <- rbind(simDat_alpha,simDat_no,simDat_beta,simDat_both)
for (s in 1:Nsub) {
  for (gen in models) {
    print(paste("Running participant",s,"/",Nsub))
    temp_dat <- subset(simDat,sub==subs[s]&gen_model==gen)
    
    ldc_file <- paste0('ldc_nn/recovery/evidence_known/sim_',gen,'_learn/no_model/ldcfit_',subs[s],'.Rdata')
    if (file.exists(ldc_file)) {
      load(ldc_file)
    }else{
      optimal_params <- DEoptim(ldc.nn.fit.w,ddm_params=ddm_params, 
                                obs = temp_dat,
                                lower = c(0,-100,1,0,0), 
                                upper = c(100,100,1,0,0),
                                Nupdate_per_trial=Nupdate_per_trial,
                                dt = dt, sigma = sigma,binning=binning,eta_sep=T,estimate_evidence=F,
                                Nsim_err=Nsim_err,error_type1=error_type1,fitname=fitname,
                                error_type2=error_type2,targetname=target,cost="separated",
                                control=c(itermax=1000,steptol=70,NP=20))
      ldc.results <- summary(optimal_params)
      #save individual results
      save(ldc.results, file=ldc_file)
    }
    fit_row <- fit_par_ev_known$sub==subs[s]&fit_par_ev_known$gen_model==gen&fit_par_ev_known$fit_model=="no"
    fit_par_ev_known[fit_row,"a0"] <- ldc.results$optim$bestmem[1]
    fit_par_ev_known[fit_row,"b0"] <- ldc.results$optim$bestmem[2]
    fit_par_ev_known[fit_row,"eta_a"] <- ldc.results$optim$bestmem[4]
    fit_par_ev_known[fit_row,"eta_b"] <- ldc.results$optim$bestmem[5]
    fit_par_ev_known[fit_row,"cost_ldc"] <- ldc.results$optim$bestval
    
    ldc_file <- paste0('ldc_nn/recovery/evidence_known/sim_',gen,'_learn/alpha_model/ldcfit_',subs[s],'.Rdata')
    if (file.exists(ldc_file)) {
      load(ldc_file)
    }else{
      optimal_params <- DEoptim(ldc.nn.fit.w,ddm_params=ddm_params, 
                                obs = temp_dat,
                                lower = c(0,-100,1,0,0), 
                                upper = c(100,100,1,10000,0),
                                Nupdate_per_trial=Nupdate_per_trial,
                                dt = dt, sigma = sigma,binning=binning,eta_sep=T,estimate_evidence=F,
                                Nsim_err=Nsim_err,error_type1=error_type1,fitname=fitname,
                                error_type2=error_type2,targetname=target,cost="separated",
                                control=c(itermax=1000,steptol=70,NP=30))
      ldc.results <- summary(optimal_params)
      #save individual results
      save(ldc.results, file=ldc_file)
    }
    fit_row <- fit_par_ev_known$sub==subs[s]&fit_par_ev_known$gen_model==gen&fit_par_ev_known$fit_model=="alpha"
    fit_par_ev_known[fit_row,"a0"] <- ldc.results$optim$bestmem[1]
    fit_par_ev_known[fit_row,"b0"] <- ldc.results$optim$bestmem[2]
    fit_par_ev_known[fit_row,"eta_a"] <- ldc.results$optim$bestmem[4]
    fit_par_ev_known[fit_row,"eta_b"] <- ldc.results$optim$bestmem[5]
    fit_par_ev_known[fit_row,"cost_ldc"] <- ldc.results$optim$bestval
    
    ldc_file <- paste0('ldc_nn/recovery/evidence_known/sim_',gen,'_learn/beta_model/ldcfit_',subs[s],'.Rdata')
    if (file.exists(ldc_file)) {
      load(ldc_file)
    }else{
      optimal_params <- DEoptim(ldc.nn.fit.w,ddm_params=ddm_params, 
                                obs = temp_dat,
                                lower = c(0,-100,1,0,0), 
                                upper = c(100,100,1,0,10000),
                                Nupdate_per_trial=Nupdate_per_trial,
                                dt = dt, sigma = sigma,binning=binning,eta_sep=T,estimate_evidence=F,
                                Nsim_err=Nsim_err,error_type1=error_type1,fitname=fitname,
                                error_type2=error_type2,targetname=target,cost="separated",
                                control=c(itermax=1000,steptol=70,NP=30))
      ldc.results <- summary(optimal_params)
      #save individual results
      save(ldc.results, file=ldc_file)
    }
    fit_row <- fit_par_ev_known$sub==subs[s]&fit_par_ev_known$gen_model==gen&fit_par_ev_known$fit_model=="beta"
    fit_par_ev_known[fit_row,"a0"] <- ldc.results$optim$bestmem[1]
    fit_par_ev_known[fit_row,"b0"] <- ldc.results$optim$bestmem[2]
    fit_par_ev_known[fit_row,"eta_a"] <- ldc.results$optim$bestmem[4]
    fit_par_ev_known[fit_row,"eta_b"] <- ldc.results$optim$bestmem[5]
    fit_par_ev_known[fit_row,"cost_ldc"] <- ldc.results$optim$bestval
    
    ldc_file <- paste0('ldc_nn/recovery/evidence_known/sim_',gen,'_learn/both_model/ldcfit_',subs[s],'.Rdata')
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
    fit_row <- fit_par_ev_known$sub==subs[s]&fit_par_ev_known$gen_model==gen&fit_par_ev_known$fit_model=="both"
    fit_par_ev_known[fit_row,"a0"] <- ldc.results$optim$bestmem[1]
    fit_par_ev_known[fit_row,"b0"] <- ldc.results$optim$bestmem[2]
    fit_par_ev_known[fit_row,"eta_a"] <- ldc.results$optim$bestmem[4]
    fit_par_ev_known[fit_row,"eta_b"] <- ldc.results$optim$bestmem[5]
    fit_par_ev_known[fit_row,"cost_ldc"] <- ldc.results$optim$bestval
  }
  fit_par_ev_known[fit_par_ev_known$sub==subs[s],"manip"] <- unique(temp_dat$manip)
}

# Model recovery ----------------------------------------------------------

###
#' List of things to do :
#' - Check that the simulated dataset were indeed generated by the intended parameters
#' - Retrieve model fits
#' - Check that each model fit is correct (i.e. correct model was called)
#' - Model comparison matrix (BIC and correlation)
#' - Is the winning model the same as the generative one ?
###

### First check that everything went fine
# Correct generative parameters
with(par,aggregate(eta_a,list(model),mean))
with(par,aggregate(eta_b,list(model),mean))
# Correct models retrieved
with(fit_par,aggregate(eta_a,list(gen=gen_model,fit=fit_model),mean))
with(fit_par,aggregate(eta_b,list(gen=gen_model,fit=fit_model),mean))
with(fit_par_ev_known,aggregate(eta_a,list(gen=gen_model,fit=fit_model),mean))
with(fit_par_ev_known,aggregate(eta_b,list(gen=gen_model,fit=fit_model),mean))
# All participants were retrieved
table(complete.cases(fit_par$a0))
table(complete.cases(fit_par_ev_known$a0))

bic_custom <- function(Residuals,k,n){
  return(log(n)*k+n*log(Residuals/n))
}

# Model comparison matrix - Evidence unknown
fit_par$bic <- bic_custom(fit_par$cost_ldc,fit_par$Npar,fit_par$Ndata_point)
mean_bic <- with(fit_par,aggregate(bic,by=list(fit=fit_model,gen=gen_model),mean))
mean_bic <- cast(mean_bic,fit~gen)

bic_sub <- with(fit_par,aggregate(bic,by=list(fit=fit_model,gen=gen_model,sub=sub),mean))
bic_sub <- cast(bic_sub,gen+sub~fit)
bic_sub$win_model <- sort(models)[apply(bic_sub[,3:6],1,which.min)]
with(bic_sub,aggregate(win_model,list(gen=gen),table))
table(subset(bic_sub,gen=="no")$win_model)
table(subset(bic_sub,gen=="alpha")$win_model)
table(subset(bic_sub,gen=="beta")$win_model)
table(subset(bic_sub,gen=="both")$win_model)

# Model comparison matrix - Evidence known
fit_par_ev_known$bic <- bic_custom(fit_par_ev_known$cost_ldc,fit_par_ev_known$Npar,fit_par_ev_known$Ndata_point)
mean_bic_ev_known <- with(fit_par_ev_known,aggregate(bic,by=list(fit=fit_model,gen=gen_model),mean))
mean_bic_ev_known <- cast(mean_bic_ev_known,fit~gen)

bic_sub_ev_known <- with(fit_par_ev_known,aggregate(bic,by=list(fit=fit_model,gen=gen_model,sub=sub),mean))
bic_sub_ev_known <- cast(bic_sub_ev_known,gen+sub~fit)
bic_sub_ev_known$win_model <- sort(models)[apply(bic_sub_ev_known[,3:6],1,which.min)]
with(bic_sub_ev_known,aggregate(win_model,list(gen=gen),table))
table(subset(bic_sub_ev_known,gen=="no")$win_model)
table(subset(bic_sub_ev_known,gen=="alpha")$win_model)
table(subset(bic_sub_ev_known,gen=="beta")$win_model)
table(subset(bic_sub_ev_known,gen=="both")$win_model)

# Probably remove ---------------------------------------------------------

cost_mat <- as.matrix(with(bic_sub,aggregate(win_model,list(gen=gen),table))[,2])
heatmap(t(cost_mat),Rowv = NA, Colv = NA, ylab = "Fitted model",xlab = "Generating model",
        labCol = c("alpha","beta","both","no"))

high_lr <- unique(subset(par,eta_a>5 | eta_b > 5)$sub)
mean_bic <- with(subset(fit_par,sub %in% high_lr),aggregate(bic,by=list(fit=fit_model,gen=gen_model),mean))
mean_bic <- cast(mean_bic,fit~gen)

test <- subset(fit_par,sub %in% high_lr & gen_model=='both')
bic <- with(mean_bic,aggregate(delta,list(model),mean))
names(bic) <- c('model','delta')

# bic_order <- c('no','alpha','beta','both','bin','bin_mid','ev_unknown')
# bic$model <- factor(bic$model, levels = bic_order)
# bic <- bic[order(bic$model),]
# barplot(bic[1:4,]$delta,names.arg = bic_order[1:4],xlab="Model (learning rate)", ylab = "BIC")
# barplot(bic[c(4,5,7),]$delta,names.arg = c("continuous","binned","ev unknown"),
#         xlab="Model (2 learning rates)", ylab = "BIC")
check <- subset(bic_sub,sub %in% high_lr)
with(check,aggregate(win_model,list(gen=gen),table))
fit_par$bic_diff_no <- fit_par$bic
for (s in subs) {
  for (gen in models) {
    fit_par[fit_par$sub==s&fit_par$gen_model==gen,'bic_diff_no'] <- 
      fit_par[fit_par$sub==s&fit_par$gen_model==gen,'bic_diff_no'] - 
      fit_par[fit_par$sub==s&fit_par$gen_model==gen&fit_par$fit_model=='no','bic']
  }
}
fit_par$max_lr <- apply(fit_par[,c('eta_a','eta_b')],1,max)
for (fit in models) {
  temp_par <- subset(fit_par,fit_model==fit)
  plot(temp_par$max_lr,temp_par$bic_diff_no,xlab = "Learning rate", ylab = "BIC diff with no learning model",
       bty='n', main = fit)
  abline(h=0)
  print(cor.test(temp_par$max_lr,temp_par$bic_diff_no))
  # abline(lm(temp_par$max_lr~temp_par$bic_diff_no)$coef)
}

