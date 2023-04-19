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
sourceCpp("ldc_train.cpp")
# Load dataset --------------------------------------------------------

train_alpha <- read.csv("train_alpha.csv")
train_beta <- read.csv("train_beta.csv")
train <- fastmerge(train_alpha,train_beta)

train$rt2 <- train$rt + train$RTconf
train$cj_pred <- -99
train$cj_lab <- train$cj
# train$cj <- train$cj/6

subs <- unique(train$sub)
conditions <- sort(unique(train$condition))
difflevels <- sort(unique(train$difflevel)) # Important to sort to match with drift order

# custom function to fit both conditions at the same time -----------------
ldc.fit <- function(params,ddm_params1,ddm_params2,obs1,obs2,dt=.001,sigma=0.1,
                    Nupdate_per_trial=1000,returnFit=T,
                    confRTname="RTconf",diffname="difflevel",respname="resp",
                    totRTname='rt2',targetname='cj',accname='cor',beta_input=.1,
                    error_type='mse',binning=F,nbin=6){
  fit1 <- ldc.nn.fit.w(params,obs1,ddm_params1,dt=dt,sigma=sigma,
               Nupdate_per_trial=Nupdate_per_trial,returnFit=returnFit,
               confRTname=confRTname,diffname=diffname,respname=respname,
               totRTname=totRTname,targetname=targetname,accname=accname,
               beta_input=beta_input,error_type=error_type,binning=binning,nbin=nbin)
  fit2 <- ldc.nn.fit.w(params,obs2,ddm_params2,dt=dt,sigma=sigma,
                        Nupdate_per_trial=Nupdate_per_trial,returnFit=returnFit,
                        confRTname=confRTname,diffname=diffname,respname=respname,
                        totRTname=totRTname,targetname=targetname,accname=accname,
                        beta_input=beta_input,error_type=error_type,binning=binning,nbin=nbin)
  if (returnFit) {
    return(fit1 + fit2)
  }else{
    return(c(fit1,fit2))
  }
}

# Fit dataset -------------------------------------------------------------
w0 <- c(10,0,1)
beta_input <- .1
Nupdate_per_trial <- 100
# bound, ter, z, vratio, drifts
params_lower <- c(0,0,0,1,0,0,0)
params_upper <- c(.2,2,0,1,.5,.5,.5)
dt <- .001; sigma <- .1
binning <- F

totlen <- length(conditions)*length(difflevels)*length(subs)
par <- data.frame(bound = NA, drift = NA, ter = NA, eta = NA,
                  cost_ddm = NA, cost_ldc = NA,
                  sub = rep(subs, each = length(conditions)*length(difflevels)),
                  condition = rep(conditions,each = length(difflevels),
                                  length.out=totlen),
                  difflevel = rep(difflevels, length.out = totlen),
                  manip=NA)
go_to("fit")
for (s in subs) {
  print(paste("Retrieving participant",s))
  #' FIT ALPHA EXPERIMENT
  temp_dat <- subset(train,sub==s)
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
  
  par[par$sub==s&par$condition=="minus","bound"] <- ddm.results$optim$bestmem[1]
  par[par$sub==s&par$condition=="minus","ter"] <- ddm.results$optim$bestmem[2]
  par[par$sub==s&par$condition=="minus","z"] <- ddm.results$optim$bestmem[3]
  par[par$sub==s&par$condition=="minus","vratio"] <- ddm.results$optim$bestmem[4]
  par[par$sub==s&par$condition=="minus"&par$difflevel==difflevels[1],"drift"] <- ddm.results$optim$bestmem[5]
  par[par$sub==s&par$condition=="minus"&par$difflevel==difflevels[2],"drift"] <- ddm.results$optim$bestmem[6]
  par[par$sub==s&par$condition=="minus"&par$difflevel==difflevels[3],"drift"] <- ddm.results$optim$bestmem[7]
  par[par$sub==s&par$condition=="minus","cost_ddm"] <- ddm.results$optim$bestval
  
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
  
  par[par$sub==s&par$condition=="plus","bound"] <- ddm.results$optim$bestmem[1]
  par[par$sub==s&par$condition=="plus","ter"] <- ddm.results$optim$bestmem[2]
  par[par$sub==s&par$condition=="plus","z"] <- ddm.results$optim$bestmem[3]
  par[par$sub==s&par$condition=="plus","vratio"] <- ddm.results$optim$bestmem[4]
  par[par$sub==s&par$condition=="plus"&par$difflevel==difflevels[1],"drift"] <- ddm.results$optim$bestmem[5]
  par[par$sub==s&par$condition=="plus"&par$difflevel==difflevels[2],"drift"] <- ddm.results$optim$bestmem[6]
  par[par$sub==s&par$condition=="plus"&par$difflevel==difflevels[3],"drift"] <- ddm.results$optim$bestmem[7]
  par[par$sub==s&par$condition=="plus","cost_ddm"] <- ddm.results$optim$bestval
  
  
  ldc_file <- paste0('w0/batch_',Nupdate_per_trial,'/ldcfit_',s,'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
  }else{
    optimal_params <- DEoptim(ldc.fit,ddm_params1=ddm_params1,ddm_params2=ddm_params2, 
                              obs1 = subset(temp_dat,condition=="minus"),
                              obs2 = subset(temp_dat,condition=="plus"),
                              lower = c(0,-100,1,0), 
                              upper = c(50,100,1,1),
                              Nupdate_per_trial=Nupdate_per_trial,
                              dt = dt, sigma = sigma,binning=binning,
                              control=c(itermax=1000,steptol=20,reltol=.001,NP=30))
    ldc.results <- summary(optimal_params)
    #save individual results
    save(ldc.results, file=ldc_file)
  }
  par[par$sub==s,"a0"] <- ldc.results$optim$bestmem[1]
  par[par$sub==s,"b0"] <- ldc.results$optim$bestmem[2]
  par[par$sub==s,"eta"] <- ldc.results$optim$bestmem[4]
  par[par$sub==s,"cost_ldc"] <- ldc.results$optim$bestval
  
  train[train$sub==s&train$condition=="minus",'cj_pred'] <-
    ldc.nn.fit.w(params=ldc.results$optim$bestmem,ddm_params = ddm_params1,
               obs=subset(temp_dat,condition=="minus"),returnFit = F,
               Nupdate_per_trial=Nupdate_per_trial, binning = binning,
               dt = dt, sigma = sigma)
  train[train$sub==s&train$condition=="plus",'cj_pred'] <-
    ldc.nn.fit.w(params=ldc.results$optim$bestmem,ddm_params = ddm_params2,
                 obs=subset(temp_dat,condition=="plus"),returnFit = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma)
  
  par[par$sub==s,"manip"] <- unique(temp_dat$manip)
}

train_alpha <- subset(train,manip=="alpha")
train_beta <- subset(train,manip=="beta")
# Compute rolling mean per subject ----------------------------------------

n <- 10 # Rolling mean window size

trial_conf_sub <- with(train,aggregate(cj,by=list(condition,trial,cor,sub),mean))
names(trial_conf_sub) <- c("condition","trial","cor","sub","cj")
trial_conf_sub <- cast(trial_conf_sub, sub + cor + trial~condition)

pred_conf_sub <- with(train,aggregate(cj_pred,by=list(condition,trial,cor,sub),mean))
names(pred_conf_sub) <- c("condition","trial","cor","sub","cj")
pred_conf_sub <- cast(pred_conf_sub, sub + cor + trial~condition)

trials <- data.frame(trial=rep(0:107,each=2),cor=c(0,1),sub=rep(subs,each=108*2))

cj_ma <- merge(trial_conf_sub,trials,all=T)
cj_pred_ma <- merge(pred_conf_sub,trials,all=T)

ma <- function(x,n,names){
  return(rollapply(x[,names], width=n, FUN=function(x) mean(x, na.rm=TRUE),partial=TRUE, align="center"))
}
for (s in subs) {
  cj_ma[cj_ma$sub==s&cj_ma$cor==0,conditions] <- ma(subset(cj_ma,sub==s&cor==0),n,conditions)
  cj_ma[cj_ma$sub==s&cj_ma$cor==1,conditions] <- ma(subset(cj_ma,sub==s&cor==1),n,conditions)
  cj_pred_ma[cj_pred_ma$sub==s&cj_pred_ma$cor==0,conditions] <- ma(subset(cj_pred_ma,sub==s&cor==0),n,conditions)
  cj_pred_ma[cj_pred_ma$sub==s&cj_pred_ma$cor==1,conditions] <- ma(subset(cj_pred_ma,sub==s&cor==1),n,conditions)
}
# Plot traces -------------------------------------------------------------
se <- function(x,na.rm=F) sd(x,na.rm=na.rm)/sqrt(length(x))

width <- 16 # Plot size expressed in cm
height <- 10

go_to("plots")
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

jpeg("Alpha_data_bin.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(3,6),
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

jpeg("Alpha_pred_bin_w0.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,type='l',col='blue',ylim=c(3,6),
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



