rm(list=ls())
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
library(myPackage) # Run devtools::install_github("pledenmat/myPackage") to install this custom package
library(Rcpp)
library(reshape)
library(lmerTest); 
library(emmeans); 
library(lattice) # qqmath 
library(timeSeries) #ColSdS
library(car)
library(zoo) # rollapply
library(colorBlindness)
sourceCpp("ldc_train.cpp")
source("ldc_nn_functions.R")

plots <- F
stat_tests <- F

# Function ----------------------------------------------------------------
max_count <- function(data){
  return(max(table(data)))
}

se <- function(x,na.rm=F) sd(x,na.rm=na.rm)/sqrt(length(x))

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
# Preprocessing -----------------------------------------------
setwd(curdir)

Data <- read.csv("alternating_fb_mod.csv")
subs <- unique(Data$sub); Nsub <- length(subs)
Ntrials <- dim(Data)[1]/Nsub

# Retrieve fits -----------------------------------------------------------
beta_input <- .1
Nupdate_per_trial <- 1
model <- "allpar"
dt <- .001; sigma <- .1
binning <- F

conditions <- sort(unique(Data$condition))
difflevels <- sort(unique(Data$difflevel)) # Important to sort to match with drift order

ntrial <- Ntrials
par_trace_bfree <- data.frame(trial=0:(ntrial-1),sub=rep(subs,each=ntrial*2),
                        condition=rep(c("plus","minus"),each=ntrial),
                        alpha=NA,beta=NA)
par_trace_afree <- data.frame(trial=0:(ntrial-1),sub=rep(subs,each=ntrial*2),
                        condition=rep(c("plus","minus"),each=ntrial),
                        alpha=NA,beta=NA)
par_trace <- data.frame(trial=0:(ntrial-1),sub=rep(subs,each=ntrial*2),
                              condition=rep(c("plus","minus"),each=ntrial),
                              alpha=NA,beta=NA)


totlen <- length(conditions)*length(difflevels)*length(subs)
par <- data.frame(bound = NA, drift = NA, ter = NA, eta = NA,
                  cost_ddm = NA, cost_ldc = NA,
                  model=rep(c('afree','bfree','full'),each=length(conditions)*length(difflevels)*Nsub),
                  sub = rep(subs, each = length(conditions)*length(difflevels)),
                  condition = rep(conditions,each = length(difflevels),
                                  length.out=totlen),
                  difflevel = rep(difflevels, length.out = totlen),
                  manip=NA)
go_to("fit")
go_to("alternating_fb")
for (s in subs) {
  print(paste("Retrieving participant",s))
  #' FIT ALPHA EXPERIMENT
  temp_dat <- subset(Data,sub==s)
  ddm_file <- paste0('ddm/ddmfit_',s,'.Rdata')
  if(file.exists(ddm_file)){
    load(ddm_file)
  }else{
    next
  }
  ddm_params <- ddm.results$optim$bestmem[c(1,5:length(ddm.results$optim$bestmem))]
  
  par[par$sub==s,"bound"] <- ddm.results$optim$bestmem[1]
  par[par$sub==s,"ter"] <- ddm.results$optim$bestmem[2]
  par[par$sub==s,"z"] <- ddm.results$optim$bestmem[3]
  par[par$sub==s,"vratio"] <- ddm.results$optim$bestmem[4]
  par[par$sub==s&par$difflevel==difflevels[1],"drift"] <- ddm.results$optim$bestmem[5]
  par[par$sub==s&par$difflevel==difflevels[2],"drift"] <- ddm.results$optim$bestmem[6]
  par[par$sub==s&par$difflevel==difflevels[3],"drift"] <- ddm.results$optim$bestmem[7]
  par[par$sub==s,"cost_ddm"] <- ddm.results$optim$bestval
  
  ldc_file <- paste0('ldc_nn/batch_',Nupdate_per_trial,'/ldcfit_',s,'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
  }else{
    next  
  }
  par[par$sub==s&par$model=='full',"a0"] <- ldc.results$optim$bestmem[1]
  par[par$sub==s&par$model=='full',"b0"] <- ldc.results$optim$bestmem[2]
  par[par$sub==s&par$model=='full',"eta"] <- ldc.results$optim$bestmem[4]
  par[par$sub==s&par$model=='full',"cost_ldc"] <- ldc.results$optim$bestval
  
  results <-
    ldc.nn.fit.w(params=c(mean(par[par$sub==s&par$model=='full',"a0"]),
                          mean(par[par$sub==s&par$model=='full',"b0"]),1,
                          mean(par[par$sub==s&par$model=='full',"eta"])),
                 ddm_params = ddm_params,
                 obs=temp_dat,returnFit = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma)
  Data[Data$sub==s,'cj_pred'] <- results$pred
  
  par_trace[par_trace$sub==s,"alpha"] <-
    c(results$trace[,1],
      rep(NA,ntrial-length(results$trace[seq(Nupdate_per_trial,nrow(results$trace),Nupdate_per_trial),1])))
  par_trace[par_trace$sub==s,"beta"] <-
    c(results$trace[,2],
      rep(NA,ntrial-length(results$trace[seq(Nupdate_per_trial,nrow(results$trace),Nupdate_per_trial),2])))
  
  
  ldc_file <- paste0('ldc_nn/eta_free_a/batch_',Nupdate_per_trial,'/ldcfit_',s,'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
  }else{
    next  
  }
  par[par$sub==s&par$model=="afree","a0"] <- ldc.results$optim$bestmem[1]
  par[par$sub==s&par$model=="afree","b0"] <- ldc.results$optim$bestmem[2]
  par[par$sub==s&par$model=="afree","eta"] <- ldc.results$optim$bestmem[4]
  par[par$sub==s&par$model=="afree","cost_ldc"] <- ldc.results$optim$bestval
  
  results <-
    ldc.nn.fit.w(params=c(mean(par[par$sub==s&par$model=="afree","a0"]),
                          mean(par[par$sub==s&par$model=="afree","b0"]),1,
                          mean(par[par$sub==s&par$model=="afree","eta"])),
                 ddm_params = ddm_params,
                 obs=temp_dat,returnFit = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma)
  Data[Data$sub==s,'cj_pred_alpha'] <- results$pred
  
  par_trace_afree[par_trace_afree$sub==s,"alpha"] <-
    c(results$trace[,1],
      rep(NA,ntrial-length(results$trace[seq(Nupdate_per_trial,nrow(results$trace),Nupdate_per_trial),1])))
  par_trace_afree[par_trace_afree$sub==s,"beta"] <-
    c(results$trace[,2],
      rep(NA,ntrial-length(results$trace[seq(Nupdate_per_trial,nrow(results$trace),Nupdate_per_trial),2])))
  
  ldc_file <- paste0('ldc_nn/eta_free_b/batch_',Nupdate_per_trial,'/ldcfit_',s,'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
  }else{
    next  
  }
  par[par$sub==s&par$model=="bfree","a0"] <- ldc.results$optim$bestmem[1]
  par[par$sub==s&par$model=="bfree","b0"] <- ldc.results$optim$bestmem[2]
  par[par$sub==s&par$model=="bfree","eta"] <- ldc.results$optim$bestmem[4]
  par[par$sub==s&par$model=="bfree","cost_ldc"] <- ldc.results$optim$bestval
  
  results <-
    ldc.nn.fit.w(params=c(mean(par[par$sub==s&par$model=="bfree","a0"]),
                          mean(par[par$sub==s&par$model=="bfree","b0"]),1,
                          mean(par[par$sub==s&par$model=="bfree","eta"])),
                 ddm_params = ddm_params,
                 obs=temp_dat,returnFit = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma)
  Data[Data$sub==s,'cj_pred_beta'] <- results$pred
  
  par_trace_bfree[par_trace_bfree$sub==s,"alpha"] <-
    c(results$trace[,1],
      rep(NA,ntrial-length(results$trace[seq(Nupdate_per_trial,nrow(results$trace),Nupdate_per_trial),1])))
  par_trace_bfree[par_trace_bfree$sub==s,"beta"] <-
    c(results$trace[,2],
      rep(NA,ntrial-length(results$trace[seq(Nupdate_per_trial,nrow(results$trace),Nupdate_per_trial),2])))
  
  
  par[par$sub==s,"manip"] <- unique(temp_dat$manip)
}
summary(par)
par <- par[complete.cases(par$bound),]

# Model comparison --------------------------------------------------------
bic_custom <- function(Residuals,k,n){
  return(log(n)*k+n*log(Residuals/n))
}

par$bic <- bic_custom(par$cost_ldc,3,Ntrials)

param_nona <- par[complete.cases(par),]
mean_bic <- with(param_nona,aggregate(bic,by=list(model=model,manip=manip),mean))
mean_resid <- with(par,aggregate(cost_ldc,by=list(model=model,manip=manip),mean))

hist(subset(par,model=='bfree')$cost_ldc,breaks=seq(0,.021,length.out=20),col = rgb(0,0,1,.2),xlim = c(0,.022),
     main = 'Model comparison', xlab = "MSE")
hist(subset(par,model=='afree')$cost_ldc,breaks=seq(0,.021,length.out=20),col = rgb(0,1,0,.2),add=T)
hist(subset(par,model=='full')$cost_ldc,breaks=seq(0,.021,length.out=20),col = rgb(1,0,0,.2),add=T)
legend('topright',legend=c('eta beta > 0','eta alpha > 0','full'),bty='n',
       fill=c(rgb(0,0,1,.2),rgb(0,1,0,.2),rgb(1,0,0,.2)))

mean_bic$delta <- -99
mean_bic[mean_bic$manip=="alpha","delta"] <- 
  mean_bic[mean_bic$manip=="alpha",]$x - 
  min(mean_bic[mean_bic$manip=="alpha",]$x)
mean_bic[mean_bic$manip=="beta","delta"] <- 
  mean_bic[mean_bic$manip=="beta",]$x -
  min(mean_bic[mean_bic$manip=="beta",]$x)

par_all <- par
# Compute rolling mean per subject ----------------------------------------

n <- 10 # Rolling mean window size

trial_conf_sub <- with(Data,aggregate(cj,by=list(trial,cor,sub),mean))
names(trial_conf_sub) <- c("trial","cor","sub","cj")

pred_conf_sub <- with(Data,aggregate(cj_pred,by=list(trial,cor,sub),mean))
names(pred_conf_sub) <- c("trial","cor","sub","cj")

trials <- data.frame(trial=rep(0:(Ntrials-1),each=2),
                     cor=c(0,1),sub=rep(subs,each=Ntrials*2))

cj_ma <- merge(trial_conf_sub,trials,all=T)
cj_pred_ma <- merge(pred_conf_sub,trials,all=T)

ma <- function(x,n,names){
  return(rollapply(x[,names], width=n, FUN=function(x) mean(x, na.rm=TRUE),partial=TRUE, align="center"))
}
for (s in subs) {
  cj_ma[cj_ma$sub==s&cj_ma$cor==0,"cj"] <- ma(subset(cj_ma,sub==s&cor==0),n,"cj")
  cj_ma[cj_ma$sub==s&cj_ma$cor==1,"cj"] <- ma(subset(cj_ma,sub==s&cor==1),n,"cj")
  cj_pred_ma[cj_pred_ma$sub==s&cj_pred_ma$cor==0,"cj"] <- ma(subset(cj_pred_ma,sub==s&cor==0),n,"cj")
  cj_pred_ma[cj_pred_ma$sub==s&cj_pred_ma$cor==1,"cj"] <- ma(subset(cj_pred_ma,sub==s&cor==1),n,"cj")
}
# Plot traces -------------------------------------------------------------
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

jpeg(filename = "traces_no_outliers.jpg",units = 'cm',width = 42,height = 30,res=300)
# layout(matrix(c(1,1,3,3,2,2,4,4,5,6,7,8),ncol=3))
layout(matrix(c(1,2,9,9,5,6,11,11,1,3,10,10,5,7,12,12,1,4,13,14,5,8,15,16),ncol=3),heights = c(.05,.05,.2,.2,.05,.05,.2,.2))
par(mar=c(0,0,0,0))
plot.new()
title(cex.main=cex.title,line=title_line,main = expression(paste(alpha,"-Manipulated Feedback")))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Empirical Data"))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Model Fits"))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Weight traces"))
plot.new()
title(cex.main=cex.title,line=title_line,main = expression(paste(beta,"-Manipulated Feedback")))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Empirical Data"))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Model Fits"))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Weight traces"))
par(mar=c(4,5,0,0)+.1)
# Plot ExpA trace --------------------------------------------------------

# Remove odd subject
RTconf_rm <- subset(Data,trial>5&RTconf>5)
RTconf_rm <- unique(RTconf_rm$sub)
Data_beta <- subset(Data_beta,!(sub %in% RTconf_rm))
Data_alpha <- subset(Data_alpha,!(sub %in% RTconf_rm))

plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)
Ntrials_phase <- Ntrials/length(unique(Data$phase))

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


# Plot ExpA predictions --------------------------------------------------------
plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)

conf_min <- with(subset(cj_pred_ma,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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

# Plot ExpB trace --------------------------------------------------------
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

# Plot ExpB predictions --------------------------------------------------------
plus_first <- unique(subset(Data_beta,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_beta,phase==0&condition=='minus')$sub)

conf_min <- with(subset(cj_pred_ma,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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

# Plot parameter traces ---------------------------------------------------
par_higheta <- subset(par,eta>10000)

plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)

# Plot trace alpha experiment
alpha_trace <- with(subset(par_trace,sub %in% Data_alpha$sub),aggregate(alpha,by=list(sub=sub,trial=trial,condition=condition),mean))
alpha_trace_minus <- cast(subset(alpha_trace,sub %in% minus_first),sub ~trial, value = "x", fun.aggregate = mean)
alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first & !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
count_plus <- count_plus[2:length(count_plus)]
count_minus <- count_minus[2:length(count_minus)]

plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
     ylim = c(-10,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(alpha_trace_plus,na.rm=T),col=VERMILLION)
polygon(c(1:ntrial,ntrial:1),c(colMeans(alpha_trace_minus,na.rm=T) + 
                                 colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(alpha_trace_minus,na.rm=T) - 
                                                                                        colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus))[ntrial:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:ntrial,ntrial:1),c(colMeans(alpha_trace_plus,na.rm=T) + 
                                 colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(alpha_trace_plus,na.rm=T) - 
                                                                                      colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus))[ntrial:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))


beta_trace <- with(subset(par_trace,sub %in% Data_alpha$sub),aggregate(beta,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(0,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(beta_trace_plus,na.rm = T),col=VERMILLION)
polygon(c(1:ntrial,ntrial:1),c(colMeans(beta_trace_minus,na.rm=T) + 
                                 colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(beta_trace_minus,na.rm=T) - 
                                                                                       colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus))[ntrial:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:ntrial,ntrial:1),c(colMeans(beta_trace_plus,na.rm=T) + 
                                 colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(beta_trace_plus,na.rm=T) - 
                                                                                     colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus))[ntrial:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))

# Plot trace Beta experiment

plus_first <- unique(subset(Data_beta,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_beta,phase==0&condition=='minus')$sub)

alpha_trace <- with(subset(par_trace,sub %in% Data_beta$sub),aggregate(alpha,by=list(sub=sub,trial=trial,condition=condition),mean))
alpha_trace_minus <- cast(subset(alpha_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
count_plus <- count_plus[2:length(count_plus)]
count_minus <- count_minus[2:length(count_minus)]

plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
     ylim=c(-10,25),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(alpha_trace_plus,na.rm=T),col=VERMILLION)
polygon(c(1:ntrial,ntrial:1),c(colMeans(alpha_trace_minus,na.rm=T) + 
                                 colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(alpha_trace_minus,na.rm=T) - 
                                                                                        colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus))[ntrial:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:ntrial,ntrial:1),c(colMeans(alpha_trace_plus,na.rm=T) + 
                                 colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(alpha_trace_plus,na.rm=T) - 
                                                                                      colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus))[ntrial:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))


beta_trace <- with(subset(par_trace,sub %in% Data_beta$sub),aggregate(beta,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(0,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(beta_trace_plus,na.rm = T),col=VERMILLION)
polygon(c(1:ntrial,ntrial:1),c(colMeans(beta_trace_minus,na.rm=T) + 
                                 colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(beta_trace_minus,na.rm=T) - 
                                                                                       colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus))[ntrial:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:ntrial,ntrial:1),c(colMeans(beta_trace_plus,na.rm=T) + 
                                 colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(beta_trace_plus,na.rm=T) - 
                                                                                     colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus))[ntrial:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))

dev.off()



