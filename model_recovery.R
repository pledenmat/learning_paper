rm(list=ls())
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
library(myPackage) # Run devtools::install_github("pledenmat/myPackage") to install this custom package
library(Rcpp)
library(DEoptim)
library(ggplot2)
library(reshape)
library(fitdistrplus)
library(zoo) # rollapply
library(colorBlindness)
sourceCpp("ldc_train.cpp")
source("ldc_nn_functions.R")
set.seed(666)

ldc <- function(a,b,data){
  return( 1 /(1 + exp( (1/sqrt(data$rt2+.00001))* (-a*data$evidence2*data$resp - b))))
}

# Load data and retrieve simulation parameters ---------------------------------------------
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
if (file.exists("par_model_recovery.csv")) {
  par <- read.csv("par_model_recovery.csv")
} else {
  for (s in 1:length(subs)) {
    print(paste("Retrieving participant",s))
    #' FIT ALPHA EXPERIMENT
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
  
  par$Ndata_point <-  round(nrow(Data)/Nsub)
  par$Npar <- 3
  par[par$model=="no",'Npar'] <- 2
  par[par$model=="both",'Npar'] <- 4
  
  bic_custom <- function(Residuals,k,n){
    return(log(n)*k+n*log(Residuals/n))
  }
  
  par$bic <- bic_custom(par$cost_ldc,par$Npar,par$Ndata_point)
  model_bic <- with(par,aggregate(bic,list(sub=sub,manip=manip,model=model),mean))
  model_bic <- cast(model_bic, sub + manip ~ model)
  table(apply(model_bic[,3:5],1,which.min))
  models_ord <- names(model_bic)[3:6]
  model_bic$best <- models_ord[apply(model_bic[,3:5],1,which.min)]
  
  alpha_best <- model_bic[model_bic$best=='alpha','sub']
  par_abest <- subset(par,sub %in% alpha_best & model == 'alpha')
  beta_best <- model_bic[model_bic$best=='beta','sub']
  par_bbest <- subset(par,sub %in% beta_best & model == 'beta')
  
  m_a <- fitdist(par_abest$eta_a/(max(par_abest$eta_a)+1), distr = 'beta')
  m_b <- fitdist(par_bbest$eta_b/(max(par_bbest$eta_b)+1), distr = 'beta')
  
  eta_a <- rbeta(Nsub,shape1 = m_a$estimate[1], shape2 = m_a$estimate[2])*(max(par_abest$eta_a)+1)+1
  eta_b <- rbeta(Nsub,shape1 = m_b$estimate[1], shape2 = m_b$estimate[2])*(max(par_bbest$eta_b)+1)+1
  par$best <- NA
  for (s in 1:Nsub) {
    par[par$sub==subs[s] ,"best"] <- model_bic[model_bic$sub==subs[s],'best']
    par_alpha_learn[par_alpha_learn$sub==subs[s] ,"best"] <- model_bic[model_bic$sub==subs[s],'best']
    par_beta_learn[par_beta_learn$sub==subs[s] ,"best"] <- model_bic[model_bic$sub==subs[s],'best']
    # par[par$sub==subs[s] & par$model %in% c("alpha","both"),"eta_a"] <- eta_a[s]
    # par[par$sub==subs[s] & par$model %in% c("beta","both"),"eta_b"] <- eta_b[s]
  }
  par <- merge(par,ddm_par,by = c("sub","condition","difflevel"))
  write.csv(par,file="par_model_recovery.csv",row.names = F)
}


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
    temp_par <- subset(par,sub==subs[s]&difflevel==level)[1,]
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

#' Finally, we compute confidence.
#' 4 simulated datasets are generated, one for each model.
#' Models differ in whether alpha/beta learning rate is fixed to 0 or not
simDat_alpha <- simDat
simDat_both <- simDat
simDat_beta <- simDat
simDat_no <- simDat
for (s in 1:Nsub) {
  print(s)
  ddm_params <- subset(par,sub==subs[s])
  ddm_params <- c(mean(ddm_params$bound),
                  with(ddm_params,aggregate(drift,list(difflevel),mean))$x)
  temp_par_both <- subset(par,sub==subs[s]&model=='both')
  temp_par_alpha <- subset(par,sub==subs[s]&model=='alpha')
  temp_par_beta <- subset(par,sub==subs[s]&model=='beta')
  temp_par_no <- subset(par,sub==subs[s]&model=='no')
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
write.csv(simDat_alpha,"simDat_alpha_model_recovery.csv",row.names = F)
write.csv(simDat_both,"simDat_both_model_recovery.csv",row.names = F)
write.csv(simDat_beta,"simDat_beta_model_recovery.csv",row.names = F)
write.csv(simDat_no,"simDat_no.csv",row.names = F)
# Retrieve fits --------------------------------------------------------------

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
      ldc_file <- paste0('ldc_nn/recovery/model_recovery/sim_',gen_model,'_learn/',fit_model,'_model/ldcfit_',subs[s],'.Rdata')
      if (file.exists(ldc_file)) {
        load(ldc_file)
        fit_par[fit_par$sub==subs[s] & fit_par$gen_model==gen_model & fit_par$fit_model==fit_model,"a0"] <- ldc.results$optim$bestmem[1]
        fit_par[fit_par$sub==subs[s] & fit_par$gen_model==gen_model & fit_par$fit_model==fit_model,"b0"] <- ldc.results$optim$bestmem[2]
        fit_par[fit_par$sub==subs[s] & fit_par$gen_model==gen_model & fit_par$fit_model==fit_model,"eta_a"] <- ldc.results$optim$bestmem[4]
        fit_par[fit_par$sub==subs[s] & fit_par$gen_model==gen_model & fit_par$fit_model==fit_model,"eta_b"] <- ldc.results$optim$bestmem[5]
        fit_par[fit_par$sub==subs[s] & fit_par$gen_model==gen_model & fit_par$fit_model==fit_model,"cost_ldc"] <- ldc.results$optim$bestval
      }
      
    }
  }
  fit_par[fit_par$sub==subs[s],"manip"] <- unique(temp_dat$manip)
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
with(par,aggregate(eta_a,list(model),summary))
with(par,aggregate(eta_b,list(model),summary))
# Correct models retrieved
with(fit_par,aggregate(eta_a,list(gen=gen_model,fit=fit_model),summary))
with(fit_par,aggregate(eta_b,list(gen=gen_model,fit=fit_model),mean))
# All participants were retrieved
table(complete.cases(fit_par$a0))

bic_custom <- function(Residuals,k,n){
  return(log(n)*k+n*log(Residuals/n))
}
aic_custom <- function(Residuals,k,n){
  return(2*k+n*log(Residuals/n))
}

fit_par$bic <- bic_custom(fit_par$cost_ldc,fit_par$Npar,fit_par$Ndata_point)
mean_bic <- with(fit_par,aggregate(bic,by=list(fit=fit_model,gen=gen_model),mean))
mean_bic <- cast(mean_bic,fit~gen)

bic_sub <- with(fit_par,aggregate(bic,by=list(fit=fit_model,gen=gen_model,sub=sub),mean))
bic_sub <- cast(bic_sub,gen+sub~fit)
bic_sub <- bic_sub[complete.cases(bic_sub$alpha),]
bic_sub$win_model <- sort(models)[apply(bic_sub[,3:6],1,which.min)]
bic_sub$worst_model <- sort(models)[apply(bic_sub[,3:6],1,which.max)]
with(bic_sub,aggregate(worst_model,list(gen=gen),table))
with(bic_sub,aggregate(win_model,list(gen=gen),table))
table(subset(bic_sub,gen=='both')$win_model)
table(subset(bic_sub,gen=='alpha')$win_model)
table(subset(bic_sub,gen=='beta')$win_model)
table(subset(bic_sub,gen=='no')$win_model)
# cost_mat <- as.matrix(with(bic_sub,aggregate(win_model,list(gen=gen),table))[,2])
# heatmap(t(cost_mat),Rowv = NA, Colv = NA, ylab = "Fitted model",xlab = "Generating model",
#         labCol = c("alpha","beta","both","no"))

plot(subset(bic_sub,gen=='beta')$beta~subset(bic_sub,gen=="beta")$alpha)
test <- subset(fit_par,sub %in% high_lr & gen_model=='both')
bic <- with(mean_bic,aggregate(delta,list(model),mean))
names(bic) <- c('model','delta')
for (m in models) {
  res <- cor(subset(bic_sub,gen==m)[,3:6],method='pearson')
  print(m)
  print(round(res,2))
}
par$max_lr <- apply(par[,c('eta_a','eta_b')],1,max)
fit_par$max_gen_lr <- NA
fit_par$bic_diff_no <- fit_par$bic
for (s in subs) {
  for (gen in models) {
    fit_par[fit_par$sub==s&fit_par$gen_model==gen,'bic_diff_no'] <- 
      fit_par[fit_par$sub==s&fit_par$gen_model==gen,'bic_diff_no'] - 
      fit_par[fit_par$sub==s&fit_par$gen_model==gen&fit_par$fit_model=='no','bic']
    fit_par[fit_par$sub==s&fit_par$gen_model==gen,'max_gen_lr'] <- mean(subset(par,sub==s&model==gen)$max_lr)
    fit_par[fit_par$sub==s&fit_par$gen_model==gen,'gen_eta_a'] <- mean(subset(par,sub==s&model==gen)$eta_a)
    fit_par[fit_par$sub==s&fit_par$gen_model==gen,'gen_eta_b'] <- mean(subset(par,sub==s&model==gen)$eta_b)
  }
}
fit_par$diff_eta_a <- abs(fit_par$eta_a - fit_par$gen_eta_a)
fit_par$diff_eta_b <- abs(fit_par$eta_b - fit_par$gen_eta_b)
for (fit in models) {
  # temp_par <- subset(fit_par,fit_model==fit)
  # plot(temp_par$max_gen_lr,temp_par$bic_diff_no,xlab = "Learning rate", ylab = "BIC diff with no learning model",
  #      bty='n', main = fit)
  # abline(h=0)
  # print(cor.test(temp_par$max_gen_lr,temp_par$bic_diff_no))
  # abline(lm(temp_par$max_lr~temp_par$bic_diff_no)$coef)
  temp_par_matched <- subset(fit_par,fit_model==fit&gen_model==fit)
  # print(cor.test(temp_par_matched$eta_a,temp_par_matched$gen_eta_a))
  print(cor.test(temp_par_matched$diff_eta_b,temp_par_matched$bic_diff_no))
}

# Compute rolling mean per subject ----------------------------------------
simDat_both_old <- read.csv("simDat_both_model_recov.csv")
simDat_alpha_old <- read.csv("simDat_alpha_model_recov.csv")
simDat_beta_old <- read.csv("simDat_beta_model_recov.csv")
n <- 25 # Rolling mean window size
n_err <- 25
Ntrials <- 756
Nskip <- 0

trial_conf_sub <- with(simDat_both,aggregate(cj,by=list(trial,cor,sub),mean))
names(trial_conf_sub) <- c("trial","cor","sub","cj")

trial_conf_sub_old <- with(simDat_both_old,aggregate(cj,by=list(trial,cor,sub),mean))
names(trial_conf_sub_old) <- c("trial","cor","sub","cj")

trial_conf_alpha <- with(simDat_alpha_old,aggregate(cj,by=list(trial,cor,sub),mean))
names(trial_conf_alpha) <- c("trial","cor","sub","cj")

trial_conf_beta <- with(simDat_beta_old,aggregate(cj,by=list(trial,cor,sub),mean))
names(trial_conf_beta) <- c("trial","cor","sub","cj")

trials <- data.frame(trial=rep((1:(Ntrials))+Nskip,each=2),
                     cor=c(0,1),sub=rep(subs,each=Ntrials*2))

cj_ma <- merge(trial_conf_sub,trials,all=T)
cj_ma_old <- merge(trial_conf_sub_old,trials,all=T)
cj_ma_alpha <- merge(trial_conf_alpha,trials,all=T)
cj_ma_beta <- merge(trial_conf_beta,trials,all=T)

ma <- function(x,n,names){
  return(rollapply(x[,names], width=n, FUN=function(x) mean(x, na.rm=TRUE),partial=TRUE, align="center"))
}
for (s in subs) {
  print(s)
  cj_ma[cj_ma$sub==s&cj_ma$cor==0,"cj"] <- ma(subset(cj_ma,sub==s&cor==0),n_err,"cj")
  cj_ma[cj_ma$sub==s&cj_ma$cor==1,"cj"] <- ma(subset(cj_ma,sub==s&cor==1),n,"cj")
  cj_ma_old[cj_ma_old$sub==s&cj_ma_old$cor==0,"cj"] <- ma(subset(cj_ma_old,sub==s&cor==0),n_err,"cj")
  cj_ma_old[cj_ma_old$sub==s&cj_ma_old$cor==1,"cj"] <- ma(subset(cj_ma_old,sub==s&cor==1),n,"cj")
  cj_ma_alpha[cj_ma_alpha$sub==s&cj_ma_alpha$cor==0,"cj"] <- ma(subset(cj_ma_alpha,sub==s&cor==0),n_err,"cj")
  cj_ma_alpha[cj_ma_alpha$sub==s&cj_ma_alpha$cor==1,"cj"] <- ma(subset(cj_ma_alpha,sub==s&cor==1),n,"cj")
  cj_ma_beta[cj_ma_beta$sub==s&cj_ma_beta$cor==0,"cj"] <- ma(subset(cj_ma_beta,sub==s&cor==0),n_err,"cj")
  cj_ma_beta[cj_ma_beta$sub==s&cj_ma_beta$cor==1,"cj"] <- ma(subset(cj_ma_beta,sub==s&cor==1),n,"cj")
}

# Plot traces old simdat to assess whether correct parameters were used ---------------------------------------------------

width <- 16 # Plot size expressed in cm
height <- 10

se <- function(x,na.rm=F) sd(x,na.rm=na.rm)/sqrt(length(x))

title_line <- -2
cex.title <- 3
cex.lab <- 3
cex.axis <- 2
cex.legend <- 2
go_to("plots")
go_to("alternating_fb")
go_to("trim")
jpeg(filename = "sim_model_recovery.jpg",units = 'cm',width = 42,height = 30,res=300)

# Plot ExpA trace
plus_first <- unique(subset(simDat_both,switch==0&condition=='plus')$sub)
minus_first <- unique(subset(simDat_both,switch==0&condition=='minus')$sub)
Ntrials_phase <- Ntrials/length(unique(simDat_both$switch))

conf_min <- with(subset(cj_ma,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_plus_se) <- c("trial","cor","cj")


xlen <- dim(conf_plus)[1]/2
conf_min_err <- subset(conf_min,cor==0)$cj
conf_min_err_se <- subset(conf_min_se,cor==0)$cj
conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj
conf_plus_err <- subset(conf_plus,cor==0)$cj
conf_plus_err_se <- subset(conf_plus_se,cor==0)$cj
conf_plus_cor <- subset(conf_plus,cor==1)$cj
conf_plus_cor_se <- subset(conf_plus_se,cor==1)$cj

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.5,.9),
     main= NULL,cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
lines(conf_min_cor,col=BLUE)
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
lines(conf_plus_err,lty = 2,col=VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))
lines(conf_plus_cor,col=VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))
legend(cex=cex.legend,"bottomleft",legend = c(expression(paste("high ",alpha," feedback first")), 
                                              expression(paste("low ",alpha," feedback first"))),
       col = c(VERMILLION,BLUE), bty = 'n', lty = c(1,1))
legend(cex=cex.legend,"bottomright",legend = c("Correct trials", "Error trials"),
       bty = 'n', lty = c(1,2))

dev.off()

# Plot ExpB trace 
plus_first <- unique(subset(Data_beta,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_beta,phase==0&condition=='minus')$sub)

conf_min <- with(subset(cj_ma,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_plus_se) <- c("trial","cor","cj")

xlen <- dim(conf_plus)[1]/2
conf_min_err <- subset(conf_min,cor==0)$cj
conf_min_err_se <- subset(conf_min_se,cor==0)$cj
conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj
conf_plus_err <- subset(conf_plus,cor==0)$cj
conf_plus_err_se <- subset(conf_plus_se,cor==0)$cj
conf_plus_cor <- subset(conf_plus,cor==1)$cj
conf_plus_cor_se <- subset(conf_plus_se,cor==1)$cj

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.5,.9),     
     main= NULL,cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
lines(conf_min_cor,col=BLUE)
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
lines(conf_plus_err,lty = 2,col=VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))
lines(conf_plus_cor,col=VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))
legend(cex=cex.legend,"bottomleft",legend = c(expression(paste("high ",beta," feedback")), expression(paste("low ",beta," feedback"))),
       col = c(VERMILLION,BLUE), bty = 'n', lty = c(1,1))
dev.off()
jpeg(filename = "sim_model_recovery_old.jpg",units = 'cm',width = 42,height = 30,res=300)

# Plot ExpA trace
plus_first <- unique(subset(simDat_alpha_old,switch==0&condition=='plus')$sub)
minus_first <- unique(subset(simDat_alpha_old,switch==0&condition=='minus')$sub)
Ntrials_phase <- Ntrials/length(unique(simDat_both_old$switch))

conf_min <- with(subset(cj_ma_old,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma_old,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma_old,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma_old,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_plus_se) <- c("trial","cor","cj")


xlen <- dim(conf_plus)[1]/2
conf_min_err <- subset(conf_min,cor==0)$cj
conf_min_err_se <- subset(conf_min_se,cor==0)$cj
conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj
conf_plus_err <- subset(conf_plus,cor==0)$cj
conf_plus_err_se <- subset(conf_plus_se,cor==0)$cj
conf_plus_cor <- subset(conf_plus,cor==1)$cj
conf_plus_cor_se <- subset(conf_plus_se,cor==1)$cj

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.5,.9),
     main= NULL,cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
lines(conf_min_cor,col=BLUE)
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
lines(conf_plus_err,lty = 2,col=VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))
lines(conf_plus_cor,col=VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))
legend(cex=cex.legend,"bottomleft",legend = c(expression(paste("high ",alpha," feedback first")), 
                                              expression(paste("low ",alpha," feedback first"))),
       col = c(VERMILLION,BLUE), bty = 'n', lty = c(1,1))
legend(cex=cex.legend,"bottomright",legend = c("Correct trials", "Error trials"),
       bty = 'n', lty = c(1,2))

dev.off()
# Plot trace for simulations from alpha model (best vs no best) --------
jpeg(filename = "sim_model_recovery_alpha_not_best.jpg",units = 'cm',width = 42,height = 30,res=300)

# Plot ExpA trace
plus_first <- unique(subset(simDat_alpha_old,switch==0&condition=='plus')$sub)
minus_first <- unique(subset(simDat_alpha_old,switch==0&condition=='minus')$sub)
Ntrials_phase <- Ntrials/length(unique(simDat_both_old$switch))

conf_min <- with(subset(cj_ma_alpha,sub %in% minus_first & !(sub %in% alpha_best)),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma_alpha,sub %in% minus_first&!(sub %in% alpha_best)),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma_alpha,sub %in% plus_first& !(sub %in% alpha_best)),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma_alpha,sub %in% plus_first& !(sub %in% alpha_best)),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_plus_se) <- c("trial","cor","cj")


xlen <- dim(conf_plus)[1]/2
conf_min_err <- subset(conf_min,cor==0)$cj
conf_min_err_se <- subset(conf_min_se,cor==0)$cj
conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj
conf_plus_err <- subset(conf_plus,cor==0)$cj
conf_plus_err_se <- subset(conf_plus_se,cor==0)$cj
conf_plus_cor <- subset(conf_plus,cor==1)$cj
conf_plus_cor_se <- subset(conf_plus_se,cor==1)$cj

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.5,.9),
     main= NULL,cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
lines(conf_min_cor,col=BLUE)
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
lines(conf_plus_err,lty = 2,col=VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))
lines(conf_plus_cor,col=VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))
legend(cex=cex.legend,"bottomleft",legend = c(expression(paste("high ",alpha," feedback first")), 
                                              expression(paste("low ",alpha," feedback first"))),
       col = c(VERMILLION,BLUE), bty = 'n', lty = c(1,1))
legend(cex=cex.legend,"bottomright",legend = c("Correct trials", "Error trials"),
       bty = 'n', lty = c(1,2))

dev.off()

# Same with alpha ---------------------------------------------------------

jpeg(filename = "sim_model_recovery_beta_not_best.jpg",units = 'cm',width = 42,height = 30,res=300)

# Plot ExpA trace
plus_first <- unique(subset(simDat_beta_old,switch==0&condition=='plus')$sub)
minus_first <- unique(subset(simDat_beta_old,switch==0&condition=='minus')$sub)
Ntrials_phase <- Ntrials/length(unique(simDat_both_old$switch))

conf_min <- with(subset(cj_ma_beta,sub %in% minus_first & !(sub %in% beta_best)),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma_beta,sub %in% minus_first&!(sub %in% beta_best)),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma_beta,sub %in% plus_first& !(sub %in% beta_best)),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma_beta,sub %in% plus_first& !(sub %in% beta_best)),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_plus_se) <- c("trial","cor","cj")


xlen <- dim(conf_plus)[1]/2
conf_min_err <- subset(conf_min,cor==0)$cj
conf_min_err_se <- subset(conf_min_se,cor==0)$cj
conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj
conf_plus_err <- subset(conf_plus,cor==0)$cj
conf_plus_err_se <- subset(conf_plus_se,cor==0)$cj
conf_plus_cor <- subset(conf_plus,cor==1)$cj
conf_plus_cor_se <- subset(conf_plus_se,cor==1)$cj

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.5,.9),
     main= NULL,cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
lines(conf_min_cor,col=BLUE)
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
lines(conf_plus_err,lty = 2,col=VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))
lines(conf_plus_cor,col=VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))
legend(cex=cex.legend,"bottomleft",legend = c(expression(paste("high ",alpha," feedback first")), 
                                              expression(paste("low ",alpha," feedback first"))),
       col = c(VERMILLION,BLUE), bty = 'n', lty = c(1,1))
legend(cex=cex.legend,"bottomright",legend = c("Correct trials", "Error trials"),
       bty = 'n', lty = c(1,2))

dev.off()

p <- ggplot(subset(par,model=='alpha'), aes(x=best, y=eta_a, fill=best)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin()
p
with(subset(par,model=='alpha'),aggregate(eta_a,list(best),summary))

# Plot learning rate according to which model was best --------------------

p <- ggplot(subset(par,model=='beta'), aes(x=best, y=eta_b, fill=best)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin()
p
with(subset(par,model=='beta'),aggregate(eta_b,list(best),summary))

