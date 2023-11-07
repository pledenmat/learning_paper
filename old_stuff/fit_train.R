rm(list=ls())
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
library(myPackage)
library(gamlss.dist) # Ex-gaussian for generating confRT distribution
library(DEoptim)
source("ldc_nn_functions.R")
library(Rcpp)
library(reshape)
library(zoo) # rollapply
library(lmerTest)
sourceCpp("ldc_train.cpp")
# Load dataset --------------------------------------------------------

train_alpha <- read.csv("train_alpha.csv")
train_beta <- read.csv("train_beta.csv")
train <- fastmerge(train_alpha,train_beta)

train$rt2 <- train$rt + train$RTconf
train$cj_pred <- -99
train$cj_pred_eta0 <- -99
train$cj_pred_per_sub <- -99
train$cj_lab <- train$cj
train$cj <- train$cj/6
binning <- F
subs <- unique(train$sub)
conditions <- sort(unique(train$condition))
difflevels <- sort(unique(train$difflevel)) # Important to sort to match with drift order

# Fit dataset -------------------------------------------------------------
w0 <- c(10,0,1)
beta_input <- .1
Nupdate_per_trial <- 1000
# bound, ter, z, vratio, drifts
params_lower <- c(0,0,0,1,0,0,0)
params_upper <- c(.2,2,0,1,.5,.5,.5)
dt <- .001; sigma <- .1

totlen <- length(conditions)*length(difflevels)*length(subs)
par <- data.frame(bound = NA, drift = NA, ter = NA, eta = NA, ntrials = NA,
                            cost_ddm = NA, cost_ldc = NA, cost_ldc_eta0 = NA,
                            sub = rep(subs, each = length(conditions)*length(difflevels)),
                            condition = rep(conditions,each = length(difflevels),
                                            length.out=totlen),
                            difflevel = rep(difflevels, length.out = totlen),
                            manip=NA)
par_per_sub <- data.frame(bound = NA, drift = NA, ter = NA, eta = NA, ntrials = NA,
                  cost_ddm = NA, cost_ldc = NA, cost_ldc_eta0 = NA,
                  sub = rep(subs, each = length(conditions)*length(difflevels)),
                  condition = rep(conditions,each = length(difflevels),
                                  length.out=totlen),
                  difflevel = rep(difflevels, length.out = totlen),
                  manip=NA)
go_to("fit")
for (s in subs) {
  for (cond in conditions) {
    print(paste("Retrieving participant",s,"condition",cond))
    #' FIT ALPHA EXPERIMENT
    temp_dat <- subset(train,sub==s&condition==cond)
    ddm_file <- paste0('ddm/ddmfit_',s,'_',cond,'.Rdata')
    if(file.exists(ddm_file)){
      load(ddm_file)
    }else{
      optimal_params <- DEoptim(ddm.fit, # function to optimize
                                lower = params_lower, upper = params_upper,obs = temp_dat,
                                dt = dt, sigma = s,ntrials=20,
                                control=c(itermax=1000,steptol=50,reltol=.001,NP=50))
      ddm.results <- summary(optimal_params)
      #save individual results
      save(ddm.results, file=ddm_file)
    }
    if (cond=="minus") {
      ddm_params1 <- ddm.results$optim$bestmem[c(1,5:length(ddm.results$optim$bestmem))]
    }else if (cond == "plus") {
      ddm_params2 <- ddm.results$optim$bestmem[c(1,5:length(ddm.results$optim$bestmem))]
    }
    ddm_params <- ddm.results$optim$bestmem[c(1,5:length(ddm.results$optim$bestmem))]
    
    ldc_file <- paste0('eta/per_cond/batch_',Nupdate_per_trial,'/ldcfit_',s,'_',cond,'.Rdata')
    if (file.exists(ldc_file)) {
      load(ldc_file)
    }else{
      optimal_params <- DEoptim(ldc.nn.fit,ddm_params=ddm_params, obs = temp_dat,
                                lower = 0, upper = 1,Nupdate_per_trial=Nupdate_per_trial,
                                dt = dt, sigma = sigma,
                                control=c(itermax=1000,steptol=20,reltol=.001,NP=10))
      ldc.results <- summary(optimal_params)
      #save individual results
      save(ldc.results, file=ldc_file)
    }
    par[par$sub==s&par$condition==cond,"bound"] <- ddm.results$optim$bestmem[1]
    par[par$sub==s&par$condition==cond,"ter"] <- ddm.results$optim$bestmem[2]
    par[par$sub==s&par$condition==cond,"z"] <- ddm.results$optim$bestmem[3]
    par[par$sub==s&par$condition==cond,"vratio"] <- ddm.results$optim$bestmem[4]
    par[par$sub==s&par$condition==cond&par$difflevel==difflevels[1],"drift"] <- ddm.results$optim$bestmem[5]
    par[par$sub==s&par$condition==cond&par$difflevel==difflevels[2],"drift"] <- ddm.results$optim$bestmem[6]
    par[par$sub==s&par$condition==cond&par$difflevel==difflevels[3],"drift"] <- ddm.results$optim$bestmem[7]
    par[par$sub==s&par$condition==cond,"eta"] <- ldc.results$optim$bestmem[1]
    par[par$sub==s&par$condition==cond,"cost_ddm"] <- ddm.results$optim$bestval
    par[par$sub==s&par$condition==cond,"cost_ldc"] <- ldc.results$optim$bestval
    par[par$sub==s&par$condition==cond,"ntrials"] <- dim(temp_dat)[1]
    par_per_sub[par_per_sub$sub==s&par_per_sub$condition==cond,"ntrials"] <- dim(temp_dat)[1]
    
    train[train$sub==s&train$condition==cond,'cj_pred'] <- 
      ldc.nn.fit(params=ldc.results$optim$bestmem[1],ddm_params = ddm_params,
                 obs=temp_dat,returnFit = F,
                 w0 = w0)
    train[train$sub==s&train$condition==cond,'cj_pred_eta0'] <- 
      ldc.nn.fit(params=0,ddm_params = ddm_params,
                 obs=temp_dat,returnFit = F,
                 w0 = w0)
    
  }
  
  par[par$sub==s,"manip"] <- unique(temp_dat$manip)
  temp_dat <- subset(train,sub==s)
  
  ldc_file <- paste0('eta/batch_',Nupdate_per_trial,'/ldcfit_',s,'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
  }else{
    print(paste("Couldn't retrieve sub",s))
  }
  par_per_sub[par_per_sub$sub==s,"eta"] <- ldc.results$optim$bestmem[4]
  par_per_sub[par_per_sub$sub==s,"cost_ldc"] <- ldc.results$optim$bestval
  
  par[par$sub==s&par$condition=="minus","cost_ldc_eta0"] <- 
    ldc.nn.fit.w(params=c(w0,0),ddm_params = ddm_params1,
                 Nupdate_per_trial = Nupdate_per_trial,
                 obs=subset(temp_dat,condition=="minus"),
                 returnFit = T)
  par[par$sub==s&par$condition=="plus","cost_ldc_eta0"] <- 
    ldc.nn.fit.w(params=c(w0,0),ddm_params = ddm_params2,
                 Nupdate_per_trial = Nupdate_per_trial,
                 obs=subset(temp_dat,condition=="plus"),
                 returnFit = T)
  par[par$sub==s&par$condition=="minus","cost_ldc_per_cond"] <- 
    ldc.nn.fit.w(params=c(w0,mean(par[par$sub==s&par$condition=="minus","eta"])),
                 ddm_params = ddm_params1,
                 Nupdate_per_trial = Nupdate_per_trial,
                 obs=subset(temp_dat,condition=="minus"),
                 returnFit = T)
  par[par$sub==s&par$condition=="plus","cost_ldc_per_cond"] <- 
    ldc.nn.fit.w(params=c(w0,mean(par[par$sub==s&par$condition=="plus","eta"])),
                 ddm_params = ddm_params2,
                 Nupdate_per_trial = Nupdate_per_trial,
                 obs=subset(temp_dat,condition=="plus"),
                 returnFit = T)
  
  par_per_sub[par_per_sub$sub==s&par_per_sub$condition=="minus","cost_ldc_per_sub"] <- 
    ldc.nn.fit.w(params=c(w0,mean(par_per_sub[par_per_sub$sub==s,"eta"])),
                 ddm_params = ddm_params1,
                 Nupdate_per_trial = Nupdate_per_trial,
                 obs=subset(temp_dat,condition=="minus"),
                 returnFit = T)
  par_per_sub[par_per_sub$sub==s&par_per_sub$condition=="plus","cost_ldc_per_sub"] <- 
    ldc.nn.fit.w(params=c(w0,mean(par_per_sub[par_per_sub$sub==s,"eta"])),
                 ddm_params = ddm_params2,
                 Nupdate_per_trial = Nupdate_per_trial,
                 obs=subset(temp_dat,condition=="plus"),
                 returnFit = T)
  
  train[train$sub==s&train$condition=="minus",'cj_pred_per_sub'] <-
    ldc.nn.fit.w(params=c(w0,mean(par_per_sub[par_per_sub$sub==s,"eta"])),
                 ddm_params = ddm_params1,
                 obs=subset(temp_dat,condition=="minus"),returnFit = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma)
  train[train$sub==s&train$condition=="plus",'cj_pred_per_sub'] <-
    ldc.nn.fit.w(params=c(w0,mean(par_per_sub[par_per_sub$sub==s,"eta"])),
                 ddm_params = ddm_params2,
                 obs=subset(temp_dat,condition=="plus"),returnFit = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma)
  
  
}

train_alpha <- subset(train,manip=="alpha")
train_beta <- subset(train,manip=="beta")
# Fit only a/b ------------------------------------------------------------
rm(obs_nn)
par_ab <- data.frame(condition=rep(c("minus","plus","both"),each=length(subs)),sub=subs,
                     alpha=NA,beta=NA)
ldc <- function(data,a,b){
  return(( 1 /(1 + exp( (1/sqrt(data$rt2+.00001))* (-a*data$evidence - b)))))
}
go_to("fit")
for (s in subs) {
  temp_dat <- subset(train,sub==s)
  if ((file.exists("obs_nn.Rdata"))) {
    load("obs_nn.Rdata")
    load("fit_ab.Rdata")
  }else{
    # Get evidence prediction -------------------------------------------------
    ddm_file1 <- paste0('ddm/ddmfit_',s,'_minus.Rdata')
    if(file.exists(ddm_file1)){
      load(ddm_file1)
    }else{
      optimal_params <- DEoptim(ddm.fit, # function to optimize
                                lower = params_lower, upper = params_upper,obs = subset(temp_dat,condition=="minus"),
                                dt = dt, sigma = s,ntrials=20,
                                control=c(itermax=1000,steptol=50,reltol=.001,NP=50))
      ddm.results <- summary(optimal_params)
      #save individual results
      save(ddm.results, file=ddm_file1)
    }
    ddm_params1 <- ddm.results$optim$bestmem[c(1,5:length(ddm.results$optim$bestmem))]
    
    
    
    ddm_file2 <- paste0('ddm/ddmfit_',s,'_plus.Rdata')
    if(file.exists(ddm_file2)){
      load(ddm_file2)
    }else{
      optimal_params <- DEoptim(ddm.fit, # function to optimize
                                lower = params_lower, upper = params_upper,obs = subset(temp_dat,condition=="plus"),
                                dt = dt, sigma = s,ntrials=20,
                                control=c(itermax=1000,steptol=50,reltol=.001,NP=50))
      ddm.results <- summary(optimal_params)
      #save individual results
      save(ddm.results, file=ddm_file2)
    }
    ddm_params2 <- ddm.results$optim$bestmem[c(1,5:length(ddm.results$optim$bestmem))]
    
    
    drift <- ddm_params1[2:length(ddm_params1)] # Make it more flexible
    bound <- ddm_params1[1]
    obs1 <- subset(temp_dat,condition=="minus")
    difficulty <- sort(unique(obs1[,"difflevel"]))
    
    obs1['evidence'] <- bound
    
    obs1_nn <- obs1[rep(seq_len(nrow(obs1)), each=Nupdate_per_trial), ]
    
    for (trial in seq(1,dim(obs1_nn)[1],Nupdate_per_trial)) {
      #' Post decision drift rate sign depends on accuracy
      if (obs1_nn[trial,"cor"] %in% c(1,'correct','cor')) {
        obs1_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] <-
          obs1_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] +
          DDM_fixed_time(v = drift[difficulty==obs1_nn[trial,"difflevel"]],
                         time=obs1_nn[trial,"RTconf"],ntrials=Nupdate_per_trial,s=sigma,dt=dt)[,1]
      }else if (obs1_nn[trial,"cor"] %in% c(-1,0,'error','err')) {
        obs1_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] <-
          obs1_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] +
          DDM_fixed_time(v = - drift[difficulty==obs1_nn[trial,"difflevel"]],
                         time=obs1_nn[trial,"RTconf"],ntrials=Nupdate_per_trial,s=sigma,dt=dt)[,1]
      }
    }
    
    drift <- ddm_params2[2:length(ddm_params2)] # Make it more flexible
    bound <- ddm_params2[1]
    
    obs2 <- subset(temp_dat,condition=="plus")
    difficulty <- sort(unique(obs2[,"difflevel"]))
    
    obs2['evidence'] <- bound
    
    obs2_nn <- obs2[rep(seq_len(nrow(obs2)), each=Nupdate_per_trial), ]
    for (trial in seq(1,dim(obs2_nn)[1],Nupdate_per_trial)) {
      #' Post decision drift rate sign depends on accuracy
      if (obs2_nn[trial,"cor"] %in% c(1,'correct','cor')) {
        obs2_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] <-
          obs2_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] +
          DDM_fixed_time(v = drift[difficulty==obs2_nn[trial,"difflevel"]],
                         time=obs2_nn[trial,"RTconf"],ntrials=Nupdate_per_trial,s=sigma,dt=dt)[,1]
      }else if (obs2_nn[trial,"cor"] %in% c(-1,0,'error','err')) {
        obs2_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] <-
          obs2_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] +
          DDM_fixed_time(v = - drift[difficulty==obs2_nn[trial,"difflevel"]],
                         time=obs2_nn[trial,"RTconf"],ntrials=Nupdate_per_trial,s=sigma,dt=dt)[,1]
      }
    }
    obs1_nn$cj_pred_per_cond <- NA
    obs2_nn$cj_pred_per_cond <- NA
    obs1_nn$cj_pred <- NA
    obs2_nn$cj_pred <- NA
    
    if (!exists("obs_nn")) {
      obs_nn <- rbind(obs1_nn,obs2_nn)
    }else{
      obs_nn <- rbind(obs_nn,obs1_nn,obs2_nn)
    }
    
    
    
    
    # Fit alpha/beta ====
    
    fit <- nls(cj~( 1 /(1 + exp( (1/sqrt(rt2+.00001))* (-a*evidence - b)))),
               data = subset(obs_nn,sub==s),
               start = list(a=15,b=0),
               trace = T)
    par_ab[par_ab$condition=="both","alpha"] <- coef(fit)["a"]
    par_ab[par_ab$condition=="both","beta"] <- coef(fit)["b"]
    
    obs_nn[obs_nn$sub==s,"cj_pred"] <- ldc(obs_nn[obs_nn$sub==s,],coef(fit)["a"],coef(fit)["b"])
    
    fit <- nls(cj~( 1 /(1 + exp( (1/sqrt(rt2+.00001))* (-a*evidence - b)))),
               data = subset(obs_nn,sub==s&condition=="minus"),
               start = list(a=15,b=0),
               trace = T)
    par_ab[par_ab$condition=="minus","alpha"] <- coef(fit)["a"]
    par_ab[par_ab$condition=="minus","beta"] <- coef(fit)["b"]
    
    obs_nn[obs_nn$sub==s&obs_nn$condition=="minus","cj_pred_per_cond"] <- ldc(
      obs_nn[obs_nn$sub==s&obs_nn$condition=="minus",],coef(fit)["a"],coef(fit)["b"])
    
    fit <- nls(cj~( 1 /(1 + exp( (1/sqrt(rt2+.00001))* (-a*evidence - b)))),
               data = subset(obs_nn,sub==s&condition=="plus"),
               start = list(a=15,b=0),
               trace = T)
    par_ab[par_ab$condition=="plus","alpha"] <- coef(fit)["a"]
    par_ab[par_ab$condition=="plus","beta"] <- coef(fit)["b"]
    
    obs_nn[obs_nn$sub==s&obs_nn$condition=="plus","cj_pred_per_cond"] <- ldc(
      obs_nn[obs_nn$sub==s&obs_nn$condition=="plus",],coef(fit)["a"],coef(fit)["b"])
    
    
    
  }
}
obs_nn$mse <- (obs_nn$cj_pred - obs_nn$cj)^2
cost <- with(obs_nn,aggregate(mse,by=list(sub=sub,trial=trial),mean))
cost <- with(cost,aggregate(x,by=list(sub=sub),mean))
cost$fit <- "eta0_ab"
cost$npar <- 2

obs_nn$mse_per_cond <- (obs_nn$cj_pred_per_cond - obs_nn$cj)^2
cost_per_cond <- with(obs_nn,aggregate(mse,by=list(sub=sub,trial=trial),mean))
cost_per_cond <- with(cost_per_cond,aggregate(x,by=list(sub=sub),mean))
cost_per_cond$fit <- "eta0_ab_per_cond"
cost_per_cond$npar <- 4

# save(par_ab,file="fit_ab.Rdata")
# save(obs_nn,file="obs_nn.Rdata")
# Model comparison --------------------------------------------------------
bic_custom <- function(Residuals,k,n){
  return(log(n)*k+n*log(Residuals/n))
}
par$cost_old <- par$cost_ldc
par$cost_ldc <- par$cost_ldc_per_cond

trials <- with(subset(par,condition=="minus"),aggregate(ntrials,by=list(sub=sub),mean)) +
  with(subset(par,condition=="plus"),aggregate(ntrials,by=list(sub=sub),mean))
trials$sub <- trials$sub/2
names(trials) <- c("sub","ntrials")

par_cost <- with(subset(par,condition=="minus"),aggregate(cost_ldc,by=list(sub=sub),mean)) +
  with(subset(par,condition=="plus"),aggregate(cost_ldc,by=list(sub=sub),mean)) 
par_cost$sub <- par_cost$sub/2
par_cost$fit <- "per_cond"
par_cost$npar <- 2


par_eta0_cost <- with(subset(par,condition=="minus"),aggregate(cost_ldc_eta0,by=list(sub=sub),mean)) +
  with(subset(par,condition=="plus"),aggregate(cost_ldc_eta0,by=list(sub=sub),mean)) 
par_eta0_cost$sub <- par_eta0_cost$sub/2
par_eta0_cost$fit <- "eta0"
par_eta0_cost$npar <- 0

par_per_sub_cost <- with(par_per_sub,aggregate(cost_ldc,by=list(sub=sub),mean))
par_per_sub_cost$fit <- "per_sub"
par_per_sub_cost$npar <- 1

par_tot <- rbind(par_cost,par_eta0_cost,par_per_sub_cost,cost,cost_per_cond)
par_tot <- merge(par_tot,trials)
names(par_tot) <- c("sub","cost","fit","npar","ntrials")

par_tot$bic <- bic_custom(par_tot$cost,par_tot$npar,par_tot$ntrials)
m <- lmer(data = par_tot, bic ~ fit + (1|sub))
anova(m)
with(par_tot,aggregate(bic,by=list(fit=fit),mean))
# par$bic <- bic_custom(par$cost_ldc,1,par$ntrials)
# par$bic_eta0 <- bic_custom(par$cost_ldc_eta0,0,par$ntrials)
# par$diff <- par$bic - par$bic_eta0
# summary(par$bic - par$bic_eta0)

# Compute rolling mean per subject ----------------------------------------

n <- 10 # Rolling mean window size

trial_conf_sub <- with(train,aggregate(cj,by=list(condition,trial,cor,sub),mean))
names(trial_conf_sub) <- c("condition","trial","cor","sub","cj")
trial_conf_sub <- cast(trial_conf_sub, sub + cor + trial~condition)

pred_conf_sub <- with(train,aggregate(cj_pred,by=list(condition,trial,cor,sub),mean))
names(pred_conf_sub) <- c("condition","trial","cor","sub","cj")
pred_conf_sub <- cast(pred_conf_sub, sub + cor + trial~condition)

pred_conf_sub_eta0 <- with(train,aggregate(cj_pred_eta0,by=list(condition,trial,cor,sub),mean))
names(pred_conf_sub_eta0) <- c("condition","trial","cor","sub","cj")
pred_conf_sub_eta0 <- cast(pred_conf_sub_eta0, sub + cor + trial~condition)

pred_conf_sub_per_sub <- with(train,aggregate(cj_pred_per_sub,by=list(condition,trial,cor,sub),mean))
names(pred_conf_sub_per_sub) <- c("condition","trial","cor","sub","cj")
pred_conf_sub_per_sub <- cast(pred_conf_sub_per_sub, sub + cor + trial~condition)

trials <- data.frame(trial=rep(0:107,each=2),cor=c(0,1),sub=rep(subs,each=108*2))

cj_ma <- merge(trial_conf_sub,trials,all=T)
cj_pred_ma <- merge(pred_conf_sub,trials,all=T)
cj_pred_ma_eta0 <- merge(pred_conf_sub_eta0,trials,all=T)
cj_pred_ma_per_sub <- merge(pred_conf_sub_per_sub,trials,all=T)

ma <- function(x,n,names){
  return(rollapply(x[,names], width=n, FUN=function(x) mean(x, na.rm=TRUE),partial=TRUE, align="center"))
}
for (s in subs) {
  cj_ma[cj_ma$sub==s&cj_ma$cor==0,conditions] <- ma(subset(cj_ma,sub==s&cor==0),n,conditions)
  cj_ma[cj_ma$sub==s&cj_ma$cor==1,conditions] <- ma(subset(cj_ma,sub==s&cor==1),n,conditions)
  cj_pred_ma[cj_pred_ma$sub==s&cj_pred_ma$cor==0,conditions] <- ma(subset(cj_pred_ma,sub==s&cor==0),n,conditions)
  cj_pred_ma[cj_pred_ma$sub==s&cj_pred_ma$cor==1,conditions] <- ma(subset(cj_pred_ma,sub==s&cor==1),n,conditions)
  cj_pred_ma_eta0[cj_pred_ma_eta0$sub==s&cj_pred_ma_eta0$cor==0,conditions] <- ma(subset(cj_pred_ma_eta0,sub==s&cor==0),n,conditions)
  cj_pred_ma_eta0[cj_pred_ma_eta0$sub==s&cj_pred_ma_eta0$cor==1,conditions] <- ma(subset(cj_pred_ma_eta0,sub==s&cor==1),n,conditions)
  cj_pred_ma_per_sub[cj_pred_ma_per_sub$sub==s&cj_pred_ma_per_sub$cor==0,conditions] <- ma(subset(cj_pred_ma_per_sub,sub==s&cor==0),n,conditions)
  cj_pred_ma_per_sub[cj_pred_ma_per_sub$sub==s&cj_pred_ma_per_sub$cor==1,conditions] <- ma(subset(cj_pred_ma_per_sub,sub==s&cor==1),n,conditions)
}

# Plot traces -------------------------------------------------------------
se <- function(x,na.rm=F) sd(x,na.rm=na.rm)/sqrt(length(x))

width <- 16 # Plot size expressed in cm
height <- 10

go_to("plots")
go_to("eta")
go_to("per_cond_wrong_err")
# Plot Exp2A trace --------------------------------------------------------
conf_min <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Alpha_data.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.9),
     main= paste("Alpha empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
dev.off()
# Plot Exp2A prediction ---------------------------------------------------
conf_min <- with(subset(cj_pred_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Alpha_pred.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.9),
     main= paste("Alpha predicted confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
dev.off()

# Plot Exp2A overlap prediction --------------------------------------------------------
conf_min <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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


pred_min <- with(subset(cj_pred_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(pred_min) <- c("trial","cor","cj")
pred_min_se <- with(subset(cj_pred_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(pred_min_se) <- c("trial","cor","cj")
pred_plus <- with(subset(cj_pred_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(pred_plus) <- c("trial","cor","cj")
pred_plus_se <- with(subset(cj_pred_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
names(pred_plus_se) <- c("trial","cor","cj")

xlen <- dim(pred_plus)[1]/2
pred_min_err <- subset(pred_min,cor==0)$cj
pred_min_err_se <- subset(pred_min_se,cor==0)$cj
pred_min_cor <- subset(pred_min,cor==1)$cj
pred_min_cor_se <- subset(pred_min_se,cor==1)$cj
pred_plus_err <- subset(pred_plus,cor==0)$cj
pred_plus_err_se <- subset(pred_plus_se,cor==0)$cj
pred_plus_cor <- subset(pred_plus,cor==1)$cj
pred_plus_cor_se <- subset(pred_plus_se,cor==1)$cj


jpeg("alpha_overlap_err.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.75),     
     main= paste("alpha empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(pred_min_err,col='lightblue')
polygon(c(1:xlen,xlen:1),c(pred_min_err + pred_min_err_se,(pred_min_err - pred_min_err_se)[xlen:1]),
        border=F,col=rgb(173, 216, 230,51,maxColorValue = 255))
lines(pred_plus_err,col='orange')
polygon(c(1:xlen,xlen:1),c(pred_plus_err + pred_plus_err_se,(pred_plus_err - pred_plus_err_se)[xlen:1]),
        border=F,col=rgb(255, 165, 0,51,maxColorValue = 255))
dev.off()


jpeg("alpha_overlap_cor.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_cor,type='l',col='blue',ylim=c(.7,.9),     
     main= paste("alpha empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(pred_min_cor,col='lightblue')
polygon(c(1:xlen,xlen:1),c(pred_min_cor + pred_min_cor_se,(pred_min_cor - pred_min_cor_se)[xlen:1]),
        border=F,col=rgb(173, 216, 230,51,maxColorValue = 255))
lines(pred_plus_cor,col='orange')
polygon(c(1:xlen,xlen:1),c(pred_plus_cor + pred_plus_cor_se,(pred_plus_cor - pred_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255, 165, 0,51,maxColorValue = 255))
dev.off()

# Plot Exp2B trace --------------------------------------------------------
conf_min <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Beta_data.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.9),     
     main= paste("Beta empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
dev.off()
# Plot Exp2B prediction ---------------------------------------------------
conf_min <- with(subset(cj_pred_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Beta_pred.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.9),     
     main= paste("Beta predicted confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
dev.off()
# Plot Exp2B overlap prediction --------------------------------------------------------
conf_min <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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


pred_min <- with(subset(cj_pred_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(pred_min) <- c("trial","cor","cj")
pred_min_se <- with(subset(cj_pred_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(pred_min_se) <- c("trial","cor","cj")
pred_plus <- with(subset(cj_pred_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(pred_plus) <- c("trial","cor","cj")
pred_plus_se <- with(subset(cj_pred_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
names(pred_plus_se) <- c("trial","cor","cj")

xlen <- dim(pred_plus)[1]/2
pred_min_err <- subset(pred_min,cor==0)$cj
pred_min_err_se <- subset(pred_min_se,cor==0)$cj
pred_min_cor <- subset(pred_min,cor==1)$cj
pred_min_cor_se <- subset(pred_min_se,cor==1)$cj
pred_plus_err <- subset(pred_plus,cor==0)$cj
pred_plus_err_se <- subset(pred_plus_se,cor==0)$cj
pred_plus_cor <- subset(pred_plus,cor==1)$cj
pred_plus_cor_se <- subset(pred_plus_se,cor==1)$cj


jpeg("Beta_overlap_err.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.75),     
     main= paste("Beta empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(pred_min_err,col='lightblue')
polygon(c(1:xlen,xlen:1),c(pred_min_err + pred_min_err_se,(pred_min_err - pred_min_err_se)[xlen:1]),
        border=F,col=rgb(173, 216, 230,51,maxColorValue = 255))
lines(pred_plus_err,col='orange')
polygon(c(1:xlen,xlen:1),c(pred_plus_err + pred_plus_err_se,(pred_plus_err - pred_plus_err_se)[xlen:1]),
        border=F,col=rgb(255, 165, 0,51,maxColorValue = 255))
dev.off()


jpeg("Beta_overlap_cor.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_cor,type='l',col='blue',ylim=c(.7,.9),     
     main= paste("Beta empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(pred_min_cor,col='lightblue')
polygon(c(1:xlen,xlen:1),c(pred_min_cor + pred_min_cor_se,(pred_min_cor - pred_min_cor_se)[xlen:1]),
        border=F,col=rgb(173, 216, 230,51,maxColorValue = 255))
lines(pred_plus_cor,col='orange')
polygon(c(1:xlen,xlen:1),c(pred_plus_cor + pred_plus_cor_se,(pred_plus_cor - pred_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255, 165, 0,51,maxColorValue = 255))
dev.off()

# Plot traces fit per sub -------------------------------------------------------------
se <- function(x,na.rm=F) sd(x,na.rm=na.rm)/sqrt(length(x))

width <- 16 # Plot size expressed in cm
height <- 10

go_to("per_sub")
# Plot Exp2A trace --------------------------------------------------------
conf_min <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Alpha_data.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.9),
     main= paste("Alpha empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
dev.off()
# Plot Exp2A prediction ---------------------------------------------------
conf_min <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Alpha_pred.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.9),
     main= paste("Alpha predicted confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
dev.off()

# Plot Exp2A overlap prediction --------------------------------------------------------
conf_min <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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


pred_min <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(pred_min) <- c("trial","cor","cj")
pred_min_se <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(pred_min_se) <- c("trial","cor","cj")
pred_plus <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(pred_plus) <- c("trial","cor","cj")
pred_plus_se <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
names(pred_plus_se) <- c("trial","cor","cj")

xlen <- dim(pred_plus)[1]/2
pred_min_err <- subset(pred_min,cor==0)$cj
pred_min_err_se <- subset(pred_min_se,cor==0)$cj
pred_min_cor <- subset(pred_min,cor==1)$cj
pred_min_cor_se <- subset(pred_min_se,cor==1)$cj
pred_plus_err <- subset(pred_plus,cor==0)$cj
pred_plus_err_se <- subset(pred_plus_se,cor==0)$cj
pred_plus_cor <- subset(pred_plus,cor==1)$cj
pred_plus_cor_se <- subset(pred_plus_se,cor==1)$cj


jpeg("alpha_overlap_err.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.75),     
     main= paste("alpha empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(pred_min_err,col='lightblue')
polygon(c(1:xlen,xlen:1),c(pred_min_err + pred_min_err_se,(pred_min_err - pred_min_err_se)[xlen:1]),
        border=F,col=rgb(173, 216, 230,51,maxColorValue = 255))
lines(pred_plus_err,col='orange')
polygon(c(1:xlen,xlen:1),c(pred_plus_err + pred_plus_err_se,(pred_plus_err - pred_plus_err_se)[xlen:1]),
        border=F,col=rgb(255, 165, 0,51,maxColorValue = 255))
dev.off()


jpeg("alpha_overlap_cor.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_cor,type='l',col='blue',ylim=c(.7,.9),     
     main= paste("alpha empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(pred_min_cor,col='lightblue')
polygon(c(1:xlen,xlen:1),c(pred_min_cor + pred_min_cor_se,(pred_min_cor - pred_min_cor_se)[xlen:1]),
        border=F,col=rgb(173, 216, 230,51,maxColorValue = 255))
lines(pred_plus_cor,col='orange')
polygon(c(1:xlen,xlen:1),c(pred_plus_cor + pred_plus_cor_se,(pred_plus_cor - pred_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255, 165, 0,51,maxColorValue = 255))
dev.off()

# Plot Exp2B trace --------------------------------------------------------
conf_min <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Beta_data.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.9),     
     main= paste("Beta empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
dev.off()
# Plot Exp2B prediction ---------------------------------------------------
conf_min <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Beta_pred.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.9),     
     main= paste("Beta predicted confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
dev.off()
# Plot Exp2B overlap prediction --------------------------------------------------------
conf_min <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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


pred_min <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(pred_min) <- c("trial","cor","cj")
pred_min_se <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(pred_min_se) <- c("trial","cor","cj")
pred_plus <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(pred_plus) <- c("trial","cor","cj")
pred_plus_se <- with(subset(cj_pred_ma_per_sub,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
names(pred_plus_se) <- c("trial","cor","cj")

xlen <- dim(pred_plus)[1]/2
pred_min_err <- subset(pred_min,cor==0)$cj
pred_min_err_se <- subset(pred_min_se,cor==0)$cj
pred_min_cor <- subset(pred_min,cor==1)$cj
pred_min_cor_se <- subset(pred_min_se,cor==1)$cj
pred_plus_err <- subset(pred_plus,cor==0)$cj
pred_plus_err_se <- subset(pred_plus_se,cor==0)$cj
pred_plus_cor <- subset(pred_plus,cor==1)$cj
pred_plus_cor_se <- subset(pred_plus_se,cor==1)$cj


jpeg("Beta_overlap_err.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.75),     
     main= paste("Beta empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(pred_min_err,col='lightblue')
polygon(c(1:xlen,xlen:1),c(pred_min_err + pred_min_err_se,(pred_min_err - pred_min_err_se)[xlen:1]),
        border=F,col=rgb(173, 216, 230,51,maxColorValue = 255))
lines(pred_plus_err,col='orange')
polygon(c(1:xlen,xlen:1),c(pred_plus_err + pred_plus_err_se,(pred_plus_err - pred_plus_err_se)[xlen:1]),
        border=F,col=rgb(255, 165, 0,51,maxColorValue = 255))
dev.off()


jpeg("Beta_overlap_cor.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_cor,type='l',col='blue',ylim=c(.7,.9),     
     main= paste("Beta empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(pred_min_cor,col='lightblue')
polygon(c(1:xlen,xlen:1),c(pred_min_cor + pred_min_cor_se,(pred_min_cor - pred_min_cor_se)[xlen:1]),
        border=F,col=rgb(173, 216, 230,51,maxColorValue = 255))
lines(pred_plus_cor,col='orange')
polygon(c(1:xlen,xlen:1),c(pred_plus_cor + pred_plus_cor_se,(pred_plus_cor - pred_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255, 165, 0,51,maxColorValue = 255))
dev.off()


# Plot traces fit ab -------------------------------------------------------------
se <- function(x,na.rm=F) sd(x,na.rm=na.rm)/sqrt(length(x))

width <- 16 # Plot size expressed in cm
height <- 10

go_to("eta0")
# Plot Exp2A trace --------------------------------------------------------
conf_min <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Alpha_data.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.9),
     main= paste("Alpha empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
dev.off()
# Plot Exp2A prediction ---------------------------------------------------
conf_min <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Alpha_pred.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.9),
     main= paste("Alpha predicted confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
dev.off()

# Plot Exp2A overlap prediction --------------------------------------------------------
conf_min <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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


pred_min <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(pred_min) <- c("trial","cor","cj")
pred_min_se <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(pred_min_se) <- c("trial","cor","cj")
pred_plus <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(pred_plus) <- c("trial","cor","cj")
pred_plus_se <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
names(pred_plus_se) <- c("trial","cor","cj")

xlen <- dim(pred_plus)[1]/2
pred_min_err <- subset(pred_min,cor==0)$cj
pred_min_err_se <- subset(pred_min_se,cor==0)$cj
pred_min_cor <- subset(pred_min,cor==1)$cj
pred_min_cor_se <- subset(pred_min_se,cor==1)$cj
pred_plus_err <- subset(pred_plus,cor==0)$cj
pred_plus_err_se <- subset(pred_plus_se,cor==0)$cj
pred_plus_cor <- subset(pred_plus,cor==1)$cj
pred_plus_cor_se <- subset(pred_plus_se,cor==1)$cj


jpeg("alpha_overlap_err.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.75),     
     main= paste("alpha empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(pred_min_err,col='lightblue')
polygon(c(1:xlen,xlen:1),c(pred_min_err + pred_min_err_se,(pred_min_err - pred_min_err_se)[xlen:1]),
        border=F,col=rgb(173, 216, 230,51,maxColorValue = 255))
lines(pred_plus_err,col='orange')
polygon(c(1:xlen,xlen:1),c(pred_plus_err + pred_plus_err_se,(pred_plus_err - pred_plus_err_se)[xlen:1]),
        border=F,col=rgb(255, 165, 0,51,maxColorValue = 255))
dev.off()


jpeg("alpha_overlap_cor.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_cor,type='l',col='blue',ylim=c(.7,.9),     
     main= paste("alpha empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(pred_min_cor,col='lightblue')
polygon(c(1:xlen,xlen:1),c(pred_min_cor + pred_min_cor_se,(pred_min_cor - pred_min_cor_se)[xlen:1]),
        border=F,col=rgb(173, 216, 230,51,maxColorValue = 255))
lines(pred_plus_cor,col='orange')
polygon(c(1:xlen,xlen:1),c(pred_plus_cor + pred_plus_cor_se,(pred_plus_cor - pred_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255, 165, 0,51,maxColorValue = 255))
dev.off()

# Plot Exp2B trace --------------------------------------------------------
conf_min <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Beta_data.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.9),     
     main= paste("Beta empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
dev.off()
# Plot Exp2B prediction ---------------------------------------------------
conf_min <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Beta_pred.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.9),     
     main= paste("Beta predicted confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
dev.off()
# Plot Exp2B overlap prediction --------------------------------------------------------
conf_min <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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


pred_min <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(pred_min) <- c("trial","cor","cj")
pred_min_se <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(pred_min_se) <- c("trial","cor","cj")
pred_plus <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(pred_plus) <- c("trial","cor","cj")
pred_plus_se <- with(subset(cj_pred_ma_eta0,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
names(pred_plus_se) <- c("trial","cor","cj")

xlen <- dim(pred_plus)[1]/2
pred_min_err <- subset(pred_min,cor==0)$cj
pred_min_err_se <- subset(pred_min_se,cor==0)$cj
pred_min_cor <- subset(pred_min,cor==1)$cj
pred_min_cor_se <- subset(pred_min_se,cor==1)$cj
pred_plus_err <- subset(pred_plus,cor==0)$cj
pred_plus_err_se <- subset(pred_plus_se,cor==0)$cj
pred_plus_cor <- subset(pred_plus,cor==1)$cj
pred_plus_cor_se <- subset(pred_plus_se,cor==1)$cj


jpeg("Beta_overlap_err.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(.5,.75),     
     main= paste("Beta empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(pred_min_err,col='lightblue')
polygon(c(1:xlen,xlen:1),c(pred_min_err + pred_min_err_se,(pred_min_err - pred_min_err_se)[xlen:1]),
        border=F,col=rgb(173, 216, 230,51,maxColorValue = 255))
lines(pred_plus_err,col='orange')
polygon(c(1:xlen,xlen:1),c(pred_plus_err + pred_plus_err_se,(pred_plus_err - pred_plus_err_se)[xlen:1]),
        border=F,col=rgb(255, 165, 0,51,maxColorValue = 255))
dev.off()


jpeg("Beta_overlap_cor.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_cor,type='l',col='blue',ylim=c(.7,.9),     
     main= paste("Beta empirical confidence, rollmean of ",n,"trials"),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(pred_min_cor,col='lightblue')
polygon(c(1:xlen,xlen:1),c(pred_min_cor + pred_min_cor_se,(pred_min_cor - pred_min_cor_se)[xlen:1]),
        border=F,col=rgb(173, 216, 230,51,maxColorValue = 255))
lines(pred_plus_cor,col='orange')
polygon(c(1:xlen,xlen:1),c(pred_plus_cor + pred_plus_cor_se,(pred_plus_cor - pred_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255, 165, 0,51,maxColorValue = 255))
dev.off()



# Check DDM fit -----------------------------------------------------------
rm(Simuls)
for (s in subs) {
  for (cond in conditions) {
    if (!exists("Simuls")) {
      temp_dat <- subset(train_alpha,sub==s&condition==cond)
      temp_par1 <- subset(par,sub==s&condition==cond)
      temp_par2 <- c(colMeans(temp_dat[,c('bound','ter')]),0,1,
                     temp_dat$drift)
      Simuls <- ddm.fit(params=temp_par2,obs=temp_dat,returnFit=F,ntrials=1)
      Simuls$sub <- s
      Simuls$condition <- cond
    }else{
      temp_dat <- subset(par,sub==s&condition==cond)
      temp_par <- c(colMeans(temp_dat[,c('bound','ter')]),0,1,
                    temp_dat$drift)
      temp_sim <- ddm.fit(params=temp_par,obs=temp_dat,returnFit=F,ntrials=1)
      temp_sim$sub <- s
      temp_sim$condition <- cond
      fastmerge(Simuls,temp_sim)
    }
  }
}
# Plot DDM predictions ----------------------------------------------------

cond_ordered <- c('minus','control','plus')
## RT
for (i in 1:N) {
  layout(matrix(c(1,2,1,3,1,4),ncol=3),heights = c(.5,3))
  par(mar=c(0,0,0,0))
  plot.new()
  text(0.5,0.5,paste("subject",subs[i]),cex=cex_title,font=2)
  par(mar=c(4,3.5,1,1))
  for (c in 1:Ncond) {
    c_data <- Data[Data$cor == 1 & Data$sub==subs[i] & Data$condition==cond_ordered[c],]
    e_data <- Data[Data$cor == 0& Data$sub==subs[i]& Data$condition==cond_ordered[c],]
    e_data$rt <- - e_data$rt
    c_simul <- Simuls[Simuls$cor==1 & Simuls$sub==subs[i]& Simuls$condition==cond_ordered[c],]
    e_simul <- Simuls[Simuls$cor==0 & Simuls$sub==subs[i]& Simuls$condition==cond_ordered[c],]
    e_simul$rt <- - e_simul$rt
    
    tempC <- hist(Data$rt[Data$cor==1&Data$sub==subs[i]],breaks=seq(0,5,.2),
                  xlim=c(0,5),prob=F,col=rgb(0,1,0,.25),border="white",
                  ylab="Frequency",xlab="Reaction times (s)",cex.lab=2,
                  cex.main=1.5, cex.axis=1.5,main=cond_ordered[c])
    tempE <- hist(Data$rt[Data$cor==0&Data$sub==subs[i]],breaks=seq(0,5,.2),
                  prob=F,add=T,col=rgb(1,0,0,.25),border='white')
    Cors <- hist(Simuls$rt[Simuls$cor==1&Simuls$sub==subs[i]],breaks=seq(0,5,.2),plot=F)
    Errs <- hist(Simuls$rt[Simuls$cor==0&Simuls$sub==subs[i]],breaks=seq(0,5,.2),plot=F)
    lines(Cors$mids,Cors$density*sum(tempC$counts)/5,type='l',col='green',lwd=3)
    lines(Errs$mids,Errs$density*sum(tempE$counts)/5,type='l',col='red',lwd=3)
    legend("topright",fill=c("white","white","green","red"),border=F,
           legend=c("Simulated corrects","Simulated errors","Empirical corrects","Empirical errors"),
           col=rep(c("Green","Red"),2),bty='n',lwd=c(1,1,-1,-1))
  }
}



# For debugging/testing the functions -------------------------------------
# 
# # a,ter,z,sigma,vratio,eta,v
# params <- c(.05,.6,0,1,.1,.08,.05)
# params <- c(a,ter,z,s,vratio,.1,v)
ntrials=10;dt=.001;sigma=.1;Nupdate_per_trial=1000;estimate_evidence <- T
confRTname="RTconf";diffname="difflevel";respname="resp";
totRTname='rt2';targetname='cj';accname='cor';beta_input=.1;error_type='mse'
obs <- temp_dat
binning=T;nbin=6
params <- ldc.results$optim$bestmem[1]

# # Old confidence trace plot ---------------------------------------------
# ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)} # Rolling mean
# 
# go_to("plots")
# 
# par(mfrow=c(1,1))
# 
# # Take confidence per trial and condition
# trial_conf_alpha <- with(train_alpha,aggregate(cj,by=list(condition,trial,cor),mean))
# names(trial_conf_alpha) <- c("condition","trial","cor","cj")
# trial_conf_alpha <- cast(trial_conf_alpha, cor + trial~condition)
# 
# trial_conf_beta <- with(train_beta,aggregate(cj,by=list(condition,trial,cor),mean))
# names(trial_conf_beta) <- c("condition","trial","cor","cj")
# trial_conf_beta <- cast(trial_conf_beta, cor + trial~condition)
# 
# # Model predictions
# pred_conf_alpha <- with(train_alpha,aggregate(cj_pred,by=list(condition,trial,cor),mean))
# names(pred_conf_alpha) <- c("condition","trial","cor","cj")
# pred_conf_alpha <- cast(pred_conf_alpha, cor + trial~condition)
# 
# pred_conf_beta <- with(train_beta,aggregate(cj_pred,by=list(condition,trial,cor),mean))
# names(pred_conf_beta) <- c("condition","trial","cor","cj")
# pred_conf_beta <- cast(pred_conf_beta, cor + trial~condition)
# 
# n <- 10 # window size
# 
# # Plot alpha experiment trajectory
# jpeg("alpha_agg.jpg",width = 12,height = 8,units = 'cm',res = 300)
# plot(ma(subset(trial_conf_alpha,cor==1)$minus,n),ylim=c(.5,.9),type = 'l', 
#      col = "blue",lty=2, main = paste('Alpha training moving window of',n,'trials'), xlab = 'Trial', ylab = 'Confidence')
# lines(ma(subset(trial_conf_alpha,cor==0)$minus,n),type = 'l', col = "blue",lty=1)
# lines(ma(subset(trial_conf_alpha,cor==1)$plus,n),type = 'l', col = "red",lty=1)
# lines(ma(subset(trial_conf_alpha,cor==0)$plus,n),type = 'l', col = "red",lty=2)
# dev.off()
# #Add model predictions
# # plot(ma(subset(pred_conf_alpha,cor==1)$minus,n),ylim=c(.65,.9),type = 'l', 
# #      col = "blue",lty=2, main = paste('Alpha prediction training moving window of',n,'trials'), xlab = 'Trial', ylab = 'Confidence')
# # lines(ma(subset(pred_conf_alpha,cor==0)$minus,n),type = 'l', col = "lightblue",lty=1)
# # lines(ma(subset(pred_conf_alpha,cor==1)$plus,n),type = 'l', col = "orange",lty=1)
# # lines(ma(subset(pred_conf_alpha,cor==0)$plus,n),type = 'l', col = "orange",lty=2)
# 
# legend('topleft',legend=c("minus","plus"),lty=c(1,1),col=c('blue','red'),bty='n')
# 
# 
# # Plot beta experiment trajectory
# jpeg("beta_agg.jpg",width = 12,height = 8,units = 'cm',res = 300)
# plot(ma(subset(trial_conf_beta,cor==1)$minus,n),ylim=c(.5,.9),type = 'l', 
#      col = "lightblue",lty=2, main = paste('Beta training moving window of',n,'trials'), xlab = 'Trial', ylab = 'Confidence')
# lines(ma(subset(trial_conf_beta,cor==0)$minus,n),type = 'l', col = "blue",lty=1)
# lines(ma(subset(trial_conf_beta,cor==1)$plus,n),type = 'l', col = "red",lty=1)
# lines(ma(subset(trial_conf_beta,cor==0)$plus,n),type = 'l', col = "red",lty=2)
# dev.off()
# # Add model predictions
# # plot(ma(subset(pred_conf_beta,cor==1)$minus,n),ylim=c(.65,.9),type = 'l', 
# #      col = "lightblue",lty=2, main = paste('Beta prediction training moving window of',n,'trials'), xlab = 'Trial', ylab = 'Confidence')
# # lines(ma(subset(pred_conf_beta,cor==0)$minus,n),type = 'l', col = "lightblue",lty=1)
# # lines(ma(subset(pred_conf_beta,cor==1)$plus,n),type = 'l', col = "orange",lty=1)
# # lines(ma(subset(pred_conf_beta,cor==0)$plus,n),type = 'l', col = "orange",lty=2)
# legend('topleft',legend=c("minus","plus"),lty=c(1,1),col=c('blue','red'),bty='n')
# 

# Trace per sub --------------------------------------------------------
go_to("trace_per_sub")

for (s in subs) {
  temp_ma <- subset(cj_ma,sub==s&cor==1)
  temp_pred_ma <- subset(cj_pred_ma,sub==s&cor==1)
  
  jpeg(paste0("trace_",s,".jpg"),width=width,height=height,units = 'cm',res=300)
  plot(temp_ma$plus,type='l',col='red',ylim=c(.5,1),
       main= paste("Confidence in correct trials, sub",s),
       xlab = "Trial", ylab = "Confidence")
  lines(temp_ma$minus,col='blue')
  lines(temp_pred_ma$plus,col="orange")
  lines(temp_pred_ma$minus,col="lightblue")
  dev.off()
}
setwd("..")

