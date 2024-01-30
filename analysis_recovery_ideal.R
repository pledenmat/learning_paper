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

simDat_no <- read.csv(file='simDat_no_model_recovery_ideal.csv')
simDat_alpha <- read.csv(file='simDat_alpha_model_recovery_ideal.csv')
simDat_beta <- read.csv(file='simDat_beta_model_recovery_ideal.csv')
simDat_both <- read.csv(file='simDat_both_model_recovery_ideal.csv')
subs <- unique(simDat_alpha$sub)
Nsub <- length(subs)
# Retrieve model fits -----------------------------------------------------


models <- c("no","alpha","beta","both")

fit_par <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                      sub = rep(subs,each=length(models)^2), manip=NA, gen_model = rep(models,length(models)),
                      fit_model = rep(models,each = length(models)),Npar = NA,
                      Ndata_point = round(nrow(simDat_alpha)/Nsub))
fit_par$Npar <- 3
fit_par[fit_par$fit_model=="no",'Npar'] <- 2
fit_par[fit_par$fit_model=="both",'Npar'] <- 4
for (s in 1:Nsub) {
  temp_dat <- subset(simDat_alpha,sub==subs[s])
  print(paste("Running participant",s,"of",Nsub))
  
  for (gen_model in models) {
    for (fit_model in models) {
      ldc_file <- paste0('fit/alternating_fb/ldc_nn/recovery/model_recovery_ideal/mean_ev/sim_',gen_model,'_learn/',fit_model,'_model/ldcfit_',subs[s],'.Rdata')
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

# fit_par <- subset(fit_par,fit_model != 'both')


with(fit_par,aggregate(eta_a,list(gen=gen_model,fit=fit_model),mean))
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
# bic_sub <- bic_sub[complete.cases(bic_sub$alpha),]
bic_sub$win_model <- sort(models)[apply(bic_sub[,3:6],1,which.min)]
bic_sub$worst_model <- sort(models)[apply(bic_sub[,3:6],1,which.max)]

# Check with aic too
fit_par$aic <- aic_custom(fit_par$cost_ldc,fit_par$Npar,fit_par$Ndata_point)
mean_aic <- with(fit_par,aggregate(aic,by=list(fit=fit_model,gen=gen_model),mean))
mean_aic <- cast(mean_aic,fit~gen)

aic_sub <- with(fit_par,aggregate(aic,by=list(fit=fit_model,gen=gen_model,sub=sub),mean))
aic_sub <- cast(aic_sub,gen+sub~fit)
aic_sub <- aic_sub[complete.cases(aic_sub$alpha),]
aic_sub$win_model <- sort(models)[apply(aic_sub[,3:6],1,which.min)]
aic_sub$worst_model <- sort(models)[apply(aic_sub[,3:6],1,which.max)]

# Rolling mean -----------------------------------------------------

n <- 25 # Rolling mean window size
n_err <- 25
Ntrials <- nrow(simDat_alpha)/Nsub
Nskip <- 0

trial_conf_sub <- with(simDat_both,aggregate(cj,by=list(trial,cor,sub),mean))
names(trial_conf_sub) <- c("trial","cor","sub","cj")

trial_conf_no <- with(simDat_no,aggregate(cj,by=list(trial,cor,sub),mean))
names(trial_conf_no) <- c("trial","cor","sub","cj")

trial_conf_alpha <- with(simDat_alpha,aggregate(cj,by=list(trial,cor,sub),mean))
names(trial_conf_alpha) <- c("trial","cor","sub","cj")

trial_conf_beta <- with(simDat_beta,aggregate(cj,by=list(trial,cor,sub),mean))
names(trial_conf_beta) <- c("trial","cor","sub","cj")

trials <- data.frame(trial=rep((1:(Ntrials))+Nskip,each=2),
                     cor=c(0,1),sub=rep(subs,each=Ntrials*2))

cj_ma <- merge(trial_conf_sub,trials,all=T)
cj_ma_no <- merge(trial_conf_no,trials,all=T)
cj_ma_alpha <- merge(trial_conf_alpha,trials,all=T)
cj_ma_beta <- merge(trial_conf_beta,trials,all=T)

ma <- function(x,n,names){
  return(rollapply(x[,names], width=n, FUN=function(x) mean(x, na.rm=TRUE),partial=TRUE, align="center"))
}
for (s in subs) {
  print(s)
  cj_ma[cj_ma$sub==s&cj_ma$cor==0,"cj"] <- ma(subset(cj_ma,sub==s&cor==0),n_err,"cj")
  cj_ma[cj_ma$sub==s&cj_ma$cor==1,"cj"] <- ma(subset(cj_ma,sub==s&cor==1),n,"cj")
  cj_ma_no[cj_ma_no$sub==s&cj_ma_no$cor==0,"cj"] <- ma(subset(cj_ma_no,sub==s&cor==0),n_err,"cj")
  cj_ma_no[cj_ma_no$sub==s&cj_ma_no$cor==1,"cj"] <- ma(subset(cj_ma_no,sub==s&cor==1),n,"cj")
  cj_ma_alpha[cj_ma_alpha$sub==s&cj_ma_alpha$cor==0,"cj"] <- ma(subset(cj_ma_alpha,sub==s&cor==0),n_err,"cj")
  cj_ma_alpha[cj_ma_alpha$sub==s&cj_ma_alpha$cor==1,"cj"] <- ma(subset(cj_ma_alpha,sub==s&cor==1),n,"cj")
  cj_ma_beta[cj_ma_beta$sub==s&cj_ma_beta$cor==0,"cj"] <- ma(subset(cj_ma_beta,sub==s&cor==0),n_err,"cj")
  cj_ma_beta[cj_ma_beta$sub==s&cj_ma_beta$cor==1,"cj"] <- ma(subset(cj_ma_beta,sub==s&cor==1),n,"cj")
}


# Plot --------------------------------------------------------------------

width <- 16 # Plot size expressed in cm
height <- 10

se <- function(x,na.rm=F) sd(x,na.rm=na.rm)/sqrt(length(x))

title_line <- -2
cex.title <- 3
cex.lab <- 3
cex.axis <- 2
cex.legend <- 2
jpeg(filename = "sim_alpha_model_recovery.jpg",units = 'cm',width = 42,height = 30,res=300)
par(mfrow=c(2,1))
# Plot ExpA trace
Ntrials_phase <- Ntrials/length(unique(simDat_both$switch))

conf_min <- with(subset(cj_ma_alpha,sub==1),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma_alpha,sub==1),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")


xlen <- dim(conf_min)[1]/2
conf_min_err <- subset(conf_min,cor==0)$cj
conf_min_err_se <- subset(conf_min_se,cor==0)$cj
conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.25,.9),
     main= "alpha feedback, alpha learn",cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(conf_min_cor,col=BLUE)
legend(cex=cex.legend,"bottomright",legend = c("Correct trials", "Error trials"),
       bty = 'n', lty = c(1,2))

conf_min <- with(subset(cj_ma_alpha,sub==2),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma_alpha,sub==2),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")


xlen <- dim(conf_min)[1]/2
conf_min_err <- subset(conf_min,cor==0)$cj
conf_min_err_se <- subset(conf_min_se,cor==0)$cj
conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.25,.9),
     main= "beta feedback, alpha learn",cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(conf_min_cor,col=BLUE)
legend(cex=cex.legend,"bottomright",legend = c("Correct trials", "Error trials"),
       bty = 'n', lty = c(1,2))


dev.off()

jpeg(filename = "sim_beta_model_recovery.jpg",units = 'cm',width = 42,height = 30,res=300)
par(mfrow=c(2,1))
# Plot ExpA trace
Ntrials_phase <- Ntrials/length(unique(simDat_both$switch))

conf_min <- with(subset(cj_ma_beta,sub==1),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma_beta,sub==1),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")


xlen <- dim(conf_min)[1]/2
conf_min_err <- subset(conf_min,cor==0)$cj
conf_min_err_se <- subset(conf_min_se,cor==0)$cj
conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.25,.9),
     main= "alpha feedback, beta learn",cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(conf_min_cor,col=BLUE)
legend(cex=cex.legend,"bottomright",legend = c("Correct trials", "Error trials"),
       bty = 'n', lty = c(1,2))

conf_min <- with(subset(cj_ma_beta,sub==2),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma_beta,sub==2),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")


xlen <- dim(conf_min)[1]/2
conf_min_err <- subset(conf_min,cor==0)$cj
conf_min_err_se <- subset(conf_min_se,cor==0)$cj
conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.25,.9),
     main= "beta feedback, beta learn",cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(conf_min_cor,col=BLUE)
legend(cex=cex.legend,"bottomright",legend = c("Correct trials", "Error trials"),
       bty = 'n', lty = c(1,2))


dev.off()

jpeg(filename = "sim_no_model_recovery.jpg",units = 'cm',width = 42,height = 30,res=300)
par(mfrow=c(2,1))
# Plot ExpA trace
Ntrials_phase <- Ntrials/length(unique(simDat_both$switch))

conf_min <- with(subset(cj_ma_no,sub==1),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma_no,sub==1),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")


xlen <- dim(conf_min)[1]/2
conf_min_err <- subset(conf_min,cor==0)$cj
conf_min_err_se <- subset(conf_min_se,cor==0)$cj
conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.25,.9),
     main= "alpha feedback, no learn",cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(conf_min_cor,col=BLUE)
legend(cex=cex.legend,"bottomright",legend = c("Correct trials", "Error trials"),
       bty = 'n', lty = c(1,2))

conf_min <- with(subset(cj_ma_no,sub==2),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma_no,sub==2),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")


xlen <- dim(conf_min)[1]/2
conf_min_err <- subset(conf_min,cor==0)$cj
conf_min_err_se <- subset(conf_min_se,cor==0)$cj
conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.25,.9),
     main= "beta feedback, no learn",cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(conf_min_cor,col=BLUE)
legend(cex=cex.legend,"bottomright",legend = c("Correct trials", "Error trials"),
       bty = 'n', lty = c(1,2))


dev.off()

jpeg(filename = "sim_both_model_recovery.jpg",units = 'cm',width = 42,height = 30,res=300)
par(mfrow=c(2,1))
# Plot ExpA trace
Ntrials_phase <- Ntrials/length(unique(simDat_both$switch))

conf_min <- with(subset(cj_ma,sub==1),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub==1),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")


xlen <- dim(conf_min)[1]/2
conf_min_err <- subset(conf_min,cor==0)$cj
conf_min_err_se <- subset(conf_min_se,cor==0)$cj
conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.25,.9),
     main= "alpha feedback,  both learn",cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(conf_min_cor,col=BLUE)
legend(cex=cex.legend,"bottomright",legend = c("Correct trials", "Error trials"),
       bty = 'n', lty = c(1,2))

conf_min <- with(subset(cj_ma,sub==2),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub==2),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")


xlen <- dim(conf_min)[1]/2
conf_min_err <- subset(conf_min,cor==0)$cj
conf_min_err_se <- subset(conf_min_se,cor==0)$cj
conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.25,.9),
     main= "beta feedback, both learn",cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(conf_min_cor,col=BLUE)
legend(cex=cex.legend,"bottomright",legend = c("Correct trials", "Error trials"),
       bty = 'n', lty = c(1,2))


dev.off()

# Plot aggreg -------------------------------------------------------------


Nphase_block <- 100 
Nphase_trial <- Ntrials_phase
simDat_alpha$withinphasetrial <- simDat_alpha$trial %% Ntrials_phase 
simDat_alpha$phase_block <-  simDat_alpha$withinphasetrial %/% (Ntrials_phase/Nphase_block)
simDat_alpha$phase_block <- as.factor(simDat_alpha$phase_block)

simDat_beta$withinphasetrial <- simDat_beta$trial %% Ntrials_phase 
simDat_beta$phase_block <-  simDat_beta$withinphasetrial %/% (Ntrials_phase/Nphase_block)
simDat_beta$phase_block <- as.factor(simDat_beta$phase_block)

simDat_both$withinphasetrial <- simDat_both$trial %% Ntrials_phase 
simDat_both$phase_block <-  simDat_both$withinphasetrial %/% (Ntrials_phase/Nphase_block)
simDat_both$phase_block <- as.factor(simDat_both$phase_block)

simDat_no$withinphasetrial <- simDat_no$trial %% Ntrials_phase 
simDat_no$phase_block <-  simDat_no$withinphasetrial %/% (Ntrials_phase/Nphase_block)
simDat_no$phase_block <- as.factor(simDat_no$phase_block)

y_range <- c(.2,.95)

count_complete <- function(dat) {
  return(sum(complete.cases(dat)))
}

conf_group <- with(simDat_no,aggregate(cj,by=list(phase_block,manip,condition,cor),mean))
names(conf_group) <- c('phase_block','manip','condition','cor','cj')

trials_phase <- data.frame(phase_block=rep(0:(Nphase_block-1),each=4),
                           cor=c(0,1),condition=rep(c("minus","plus"),each=2),manip=rep(c("alpha","beta"),each=Nphase_block*4))

conf_group <- merge(conf_group,trials_phase,all=T)

conf_alpha <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)

jpeg(filename="no_sim_cj_aggreg_alpha_manip.jpg",height = 20,width=30,units = 'cm',res=600)
plot(conf_alpha$minus_0,ylim=y_range,col = BLUE, type = 'b',main="Alpha-Manipulated Feedback",
     lty = 2, pch = 17, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",cex.main=cex.title,
     xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
lines(conf_alpha$plus_0, type = 'b', pch = 17, col = VERMILLION, lwd = 2, lty = 2)
lines(conf_alpha$plus_1, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 1)
lines(conf_alpha$minus_1, type = 'b', pch = 16, col = BLUE, lwd = 2, lty = 1)
legend("top",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)
dev.off()

jpeg(filename="no_sim_cj_aggreg_beta_manip.jpg",height = 20,width=30,units = 'cm',res=600)
plot(conf_beta$minus_0,ylim=y_range,col = BLUE, type = 'b',main="Beta-Manipulated Feedback",
     lty = 2, pch = 17, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",cex.main=cex.title,
     xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
lines(conf_beta$plus_0, type = 'b', pch = 17, col = VERMILLION, lwd = 2, lty = 2)
lines(conf_beta$plus_1, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 1)
lines(conf_beta$minus_1, type = 'b', pch = 16, col = BLUE, lwd = 2, lty = 1)
legend("top",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)
dev.off()

conf_group <- with(simDat_beta,aggregate(cj,by=list(phase_block,manip,condition,cor),mean))
names(conf_group) <- c('phase_block','manip','condition','cor','cj')

trials_phase <- data.frame(phase_block=rep(0:(Nphase_block-1),each=4),
                           cor=c(0,1),condition=rep(c("minus","plus"),each=2),manip=rep(c("alpha","beta"),each=Nphase_block*4))

conf_group <- merge(conf_group,trials_phase,all=T)

conf_alpha_count <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor,count_complete)
conf_alpha <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_count <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,count_complete)

jpeg(filename="beta_sim_cj_aggreg_alpha_manip.jpg",height = 20,width=30,units = 'cm',res=600)
plot(conf_alpha$minus_0,ylim=y_range,col = BLUE, type = 'b',main="Alpha-Manipulated Feedback",
     lty = 2, pch = 17, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",cex.main=cex.title,
     xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
lines(conf_alpha$plus_0, type = 'b', pch = 17, col = VERMILLION, lwd = 2, lty = 2)
lines(conf_alpha$plus_1, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 1)
lines(conf_alpha$minus_1, type = 'b', pch = 16, col = BLUE, lwd = 2, lty = 1)
legend("top",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)
dev.off()

jpeg(filename="beta_sim_cj_aggreg_beta_manip.jpg",height = 20,width=30,units = 'cm',res=600)
plot(conf_beta$minus_0,ylim=y_range,col = BLUE, type = 'b',main="Beta-Manipulated Feedback",
     lty = 2, pch = 17, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",cex.main=cex.title,
     xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
lines(conf_beta$plus_0, type = 'b', pch = 17, col = VERMILLION, lwd = 2, lty = 2)
lines(conf_beta$plus_1, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 1)
lines(conf_beta$minus_1, type = 'b', pch = 16, col = BLUE, lwd = 2, lty = 1)
legend("top",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)
dev.off()


# Aggregating behavior
conf_group <- with(simDat_alpha,aggregate(cj,by=list(phase_block,manip,condition,cor),mean))
names(conf_group) <- c('phase_block','manip','condition','cor','cj')

trials_phase <- data.frame(phase_block=rep(0:(Nphase_block-1),each=4),
                           cor=c(0,1),condition=rep(c("minus","plus"),each=2),manip=rep(c("alpha","beta"),each=Nphase_block*4))

conf_group <- merge(conf_group,trials_phase,all=T)

conf_alpha_count <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor,count_complete)
conf_alpha <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_count <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,count_complete)

jpeg(filename="alpha_sim_cj_aggreg_alpha_manip.jpg",height = 20,width=30,units = 'cm',res=600)
plot(conf_alpha$minus_0,ylim=y_range,col = BLUE, type = 'b',main="Alpha-Manipulated Feedback",
     lty = 2, pch = 17, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",cex.main=cex.title,
     xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
lines(conf_alpha$plus_0, type = 'b', pch = 17, col = VERMILLION, lwd = 2, lty = 2)
lines(conf_alpha$plus_1, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 1)
lines(conf_alpha$minus_1, type = 'b', pch = 16, col = BLUE, lwd = 2, lty = 1)
legend("top",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)
dev.off()

jpeg(filename="alpha_sim_cj_aggreg_beta_manip.jpg",height = 20,width=30,units = 'cm',res=600)
plot(conf_beta$minus_0,ylim=y_range,col = BLUE, type = 'b',main="Beta-Manipulated Feedback",
     lty = 2, pch = 17, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",cex.main=cex.title,
     xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
lines(conf_beta$plus_0, type = 'b', pch = 17, col = VERMILLION, lwd = 2, lty = 2)
lines(conf_beta$plus_1, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 1)
lines(conf_beta$minus_1, type = 'b', pch = 16, col = BLUE, lwd = 2, lty = 1)
legend("top",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)
dev.off()


conf_group <- with(simDat_both,aggregate(cj,by=list(phase_block,manip,condition,cor),mean))
names(conf_group) <- c('phase_block','manip','condition','cor','cj')

trials_phase <- data.frame(phase_block=rep(0:(Nphase_block-1),each=4),
                           cor=c(0,1),condition=rep(c("minus","plus"),each=2),manip=rep(c("alpha","beta"),each=Nphase_block*4))

conf_group <- merge(conf_group,trials_phase,all=T)

conf_alpha_count <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor,count_complete)
conf_alpha <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_count <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,count_complete)

jpeg(filename="both_sim_cj_aggreg_alpha_manip.jpg",height = 20,width=30,units = 'cm',res=600)
plot(conf_alpha$minus_0,ylim=y_range,col = BLUE, type = 'b',main="Alpha-Manipulated Feedback",
     lty = 2, pch = 17, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",cex.main=cex.title,
     xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
lines(conf_alpha$plus_0, type = 'b', pch = 17, col = VERMILLION, lwd = 2, lty = 2)
lines(conf_alpha$plus_1, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 1)
lines(conf_alpha$minus_1, type = 'b', pch = 16, col = BLUE, lwd = 2, lty = 1)
legend("top",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)
dev.off()

jpeg(filename="both_sim_cj_aggreg_beta_manip.jpg",height = 20,width=30,units = 'cm',res=600)
plot(conf_beta$minus_0,ylim=y_range,col = BLUE, type = 'b',main="Beta-Manipulated Feedback",
     lty = 2, pch = 17, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",cex.main=cex.title,
     xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
lines(conf_beta$plus_0, type = 'b', pch = 17, col = VERMILLION, lwd = 2, lty = 2)
lines(conf_beta$plus_1, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 1)
lines(conf_beta$minus_1, type = 'b', pch = 16, col = BLUE, lwd = 2, lty = 1)
legend("top",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)
dev.off()

