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
library(emmeans)

sourceCpp("ldc_train.cpp")
# Load dataset --------------------------------------------------------

train_alpha <- read.csv("train_alpha.csv")
train_beta <- read.csv("train_beta.csv")
train <- fastmerge(train_alpha,train_beta)

train$rt2 <- train$rt + train$RTconf
train$cj_pred <- -99
train$cj_lab <- train$cj
binning <- T
train$cj <- train$cj/6
subs <- unique(train$sub)
conditions <- sort(unique(train$condition))
difflevels <- sort(unique(train$difflevel)) # Important to sort to match with drift order

# custom function to fit both conditions at the same time -----------------
ldc.fit <- function(params,ddm_params1,ddm_params2,obs1,obs2,dt=.001,sigma=0.1,
                    Nupdate_per_trial=1000,returnFit=T,
                    confRTname="RTconf",diffname="difflevel",respname="resp",
                    totRTname='rt2',targetname='cj',accname='cor',beta_input=.1,
                    error_type='mse',binning=F,nbin=6,shuffle=F){
  fit1 <- ldc.nn.fit.w(params,obs1,ddm_params1,dt=dt,sigma=sigma,
                       Nupdate_per_trial=Nupdate_per_trial,returnFit=returnFit,
                       confRTname=confRTname,diffname=diffname,respname=respname,
                       totRTname=totRTname,targetname=targetname,accname=accname,
                       beta_input=beta_input,error_type=error_type,binning=binning,nbin=nbin,shuffle=shuffle)
  fit2 <- ldc.nn.fit.w(params,obs2,ddm_params2,dt=dt,sigma=sigma,
                       Nupdate_per_trial=Nupdate_per_trial,returnFit=returnFit,
                       confRTname=confRTname,diffname=diffname,respname=respname,
                       totRTname=totRTname,targetname=targetname,accname=accname,
                       beta_input=beta_input,error_type=error_type,binning=binning,nbin=nbin,shuffle=shuffle)
  if (returnFit) {
    return(fit1 + fit2)
  }else{
    return(c(fit1,fit2))
  }
}

# Fit dataset -------------------------------------------------------------
w0 <- c(10,0,1)
beta_input <- .1
Nupdate_per_trial <- 500
# bound, ter, z, vratio, drifts
params_lower <- c(0,0,0,1,0,0,0)
params_upper <- c(.2,2,0,1,.5,.5,.5)
dt <- .001; sigma <- .1

totlen <- length(conditions)*length(difflevels)*length(subs)
par <- data.frame(eta = NA, a0 = NA, b0 = NA, cost_ldc = NA,
                  sub = rep(subs, each = length(conditions)*length(difflevels)),
                  condition = rep(conditions,each = length(difflevels),
                                  length.out=totlen),
                  difflevel = rep(difflevels, length.out = totlen),
                  manip=NA)
par2 <- data.frame(eta = NA, a0 = NA, b0 = NA, cost_ldc = NA,
                   sub = rep(subs, each = length(conditions)*length(difflevels)),
                   condition = rep(conditions,each = length(difflevels),
                                   length.out=totlen),
                   difflevel = rep(difflevels, length.out = totlen),
                   manip=NA)
go_to("fit")
for (s in subs) {
  print(paste("Retrieving participant",s))
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
  
  
  # ldc_file <- paste0('error_integrate/batch_',Nupdate_per_trial,'/ldcfit_',s,'.Rdata')
  # if (file.exists(ldc_file)) {
  #   load(ldc_file)
  # }else{
  #   print(paste("Failed to retrieve participant",s))
  #   next
  # }
  # par[par$sub==s,"a0"] <- ldc.results$optim$bestmem[1]
  # par[par$sub==s,"b0"] <- ldc.results$optim$bestmem[2]
  # par[par$sub==s,"eta"] <- ldc.results$optim$bestmem[4]
  # par[par$sub==s,"cost_ldc"] <- ldc.results$optim$bestval
  # par[par$sub==s,"manip"] <- unique(temp_dat$manip)
  
  ldc_file <- paste0('error_integrate/batch_',Nupdate_per_trial,'/ldcfit_',s,'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
  }else{
    print(paste("Failed to retrieve participant",s))
  }
  par2[par2$sub==s,"a0"] <- ldc.results$optim$bestmem[1]
  par2[par2$sub==s,"b0"] <- ldc.results$optim$bestmem[2]
  par2[par2$sub==s,"eta"] <- ldc.results$optim$bestmem[4]
  par2[par2$sub==s,"cost_ldc"] <- ldc.results$optim$bestval
  par2[par2$sub==s,"manip"] <- unique(temp_dat$manip)
  
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
}

train_alpha <- subset(train,manip=="alpha")
train_beta <- subset(train,manip=="beta")


# par$fit <- 1000
# par2$fit <- 500

# par <- par[complete.cases(par),]
# par2 <- par2[complete.cases(par2),]
# par_tot <- rbind(par,par2)

# m <- lmer(data = par_tot, eta ~ fit + (1|sub))
# anova(m)
# emm <- emmeans(m, ~ fit)
# pairs(emm)
# m <- lmer(data = par_tot, a0 ~ fit + (1|sub))
# anova(m)
# emm <- emmeans(m, ~ fit)
# pairs(emm)
# m <- lmer(data = par_tot, b0 ~ fit + (1|sub))
# anova(m)
# emm <- emmeans(m, ~ fit)
# pairs(emm)
# with(par_tot,aggregate(a0,by=list(fit),mean))
# with(par_tot,aggregate(b0,by=list(fit),mean))
# with(par_tot,aggregate(eta,by=list(fit),mean))
# 
# cor.test(par$eta,par2$eta)
# cor.test(par$a0,par2$a0)
# cor.test(par$b0,par2$b0)
# cor.test(par$a0,parfast_shuffle$a0)
# 
# m <- lmer(data = par_tot, cost_ldc ~ fit + (1|sub))
# anova(m)
# emm <- emmeans(m, ~ fit)
# pairs(emm)
# with(par_tot,aggregate(cost_ldc,by=list(fit),mean))
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
go_to("memory")
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
plot(conf_min_err,bty='n',lty = 2,type='l',col='blue',ylim=c(.5,.9),
     main= expression(paste("Exp ",alpha," empirical confidence")),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,lty = 2,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
legend("bottomleft",legend = c(expression(paste("high ",alpha," feedback")), expression(paste("low ",alpha," feedback"))),
       col = c("red","blue"), bty = 'n', lty = c(1,1))
legend("bottomright",legend = c("Correct trials", "Error trials"),
       bty = 'n', lty = c(1,2))
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
plot(conf_min_err,bty='n',lty = 2,type='l',col='blue',ylim=c(.5,.9),
     main= expression(paste("Exp ",alpha," predicted confidence")),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,lty = 2,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
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
plot(conf_min_err,bty='n',lty = 2,type='l',col='blue',ylim=c(.5,.9),     
     main= expression(paste("Exp ",beta," empirical confidence")),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,lty = 2,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
legend("bottomleft",legend = c(expression(paste("high ",beta," feedback")), expression(paste("low ",beta," feedback"))),
       col = c("red","blue"), bty = 'n', lty = c(1,1))
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
plot(conf_min_err,bty='n',lty = 2,type='l',col='blue',ylim=c(.5,.9),     
     main= expression(paste("Exp ",beta," predicted confidence")),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,lty = 2,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
dev.off()
