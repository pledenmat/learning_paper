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

# Setup simulation parameters ---------------------------------------------
### DDM parameters
bound <- .05
drift <- c(.03,.05,.07)
difflevels <- c('hard','medium','easy')
ter <- .4
conf_ter <- .4
z <- 0
sigma <- .1
dt <- .001

### LDC parameters
a0 <- 15
b0 <- 0
eta_a <- c(0,5,0,5)
eta_b <- c(0,0,5,5)
models <- c('no','alpha','beta','both')
Nupdate_per_trial <- 1
binning <- F

### Experiment
fb_manip <- c('alpha','beta')
alpha_minus <- 1
alpha_neutral <- 18
alpha_plus <- 50
beta_minus <- -2
beta_neutral <- 0
beta_plus <- 2
nswitch <- 100
ntrial_per_switch <- 100 # Actually multiplied by number of difficulty levels (3)
ntrial <- nswitch*ntrial_per_switch

#' First, we generate trials from DDM parameters
rm(simDat)
for (manip in 1:length(fb_manip)) {
  for(diff in 1:length(difflevels)){
    predictions <- data.frame(DDM_with_confidence_slow(
      v=drift[diff],a=bound,ter=ter,z=z,ntrials=ntrial,
      s=sigma,dt=dt,t2time = 0,postdriftmod = 1))
    names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
    predictions$RTconf <- predictions$rt - ter + conf_ter
    predictions$RTconf <- sample(predictions$RTconf)
    predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(
      v=drift[diff],a=bound,ter=ter,z=z,ntrials=ntrial,
      s=sigma,dt=dt,t2distribution = predictions$RTconf,postdriftmod = 1))
    names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
    predictions$RTconf <- predictions$rt2 - predictions$rt
    predictions$difflevel <- difflevels[diff]
    predictions$manip <- fb_manip[manip]
    predictions$sub <- manip
    predictions$switch <- 1:nrow(predictions) %% nswitch
    predictions$trial <- -99 # Actual trial numbers are given at the subject level below
    predictions$condition <- "minus"
    predictions[predictions$switch %% 2 == 1,"condition"] <- "plus"
    if(diff==1&manip==1){ simDat <- predictions
    }else{
      simDat <- rbind(simDat,predictions)
    }
  }
  #' We reorganize trial order to reproduce the actual experiment, i.e.:
  #' - Intermixed trial difficulty
  #' - Alternating feedback conditions
  temp <- simDat[simDat$sub==manip,]
  temp <- temp[order(temp$sub,temp$switch),]
  count_switch <- table(temp$switch)
  cum_switch <- cumsum(count_switch)
  shuffle <- sample(1:cum_switch[1])
  for (iter_switch in 1:(nswitch-1)) {
    shuffle <- c(shuffle,sample((cum_switch[iter_switch]+1):cum_switch[iter_switch+1]))
  }
  simDat[simDat$sub==manip,] <- temp[shuffle,]
  simDat[simDat$sub==manip,"trial"] <- 1:nrow(simDat[simDat$sub==manip,])
}


# Then we provide feedback from same parameters than the experiment
simDat$fb <- -99
simDat[simDat$manip=='alpha'&simDat$condition=='minus','fb'] <- ldc(alpha_minus,beta_neutral,simDat[simDat$manip=='alpha'&simDat$condition=='minus',])
simDat[simDat$manip=='beta'&simDat$condition=='minus','fb'] <- ldc(alpha_neutral,beta_minus,simDat[simDat$manip=='beta'&simDat$condition=='minus',])
simDat[simDat$manip=='alpha'&simDat$condition=='plus','fb'] <- ldc(alpha_plus,beta_neutral,simDat[simDat$manip=='alpha'&simDat$condition=='plus',])
simDat[simDat$manip=='beta'&simDat$condition=='plus','fb'] <- ldc(alpha_neutral,beta_plus,simDat[simDat$manip=='beta'&simDat$condition=='plus',])

simDat$evidence <- simDat$evidence2 # Variable naming for generating confidence

simDat$evidence_mean_ev <- bound # Variable naming for generating confidence
# Add a drift rate column
simDat$drift <- drift[1]
for (diff in 2:length(drift)) {
  simDat[simDat$difflevel==difflevels[diff],'drift'] <- drift[diff]
}
# Add evidence if correct, deduce if error
simDat$evidence_mean_ev <- simDat$evidence_mean_ev + (as.numeric(simDat$cor)-.5)*2 * simDat$drift * simDat$RTconf 

#' Finally, we compute confidence.
#' 4 simulated datasets are generated, one for each model.
#' Models differ in whether alpha/beta learning rate is fixed to 0 or not
simDat_alpha <- simDat
simDat_both <- simDat
simDat_beta <- simDat
simDat_no <- simDat
for (iter_manip in 1:length(fb_manip)) {
  temp_dat <- subset(simDat,sub==iter_manip)
  results_alpha_learn <-
    ldc.nn.fit.w(params=c(a0,b0,1,eta_a[which(models=='alpha')],
                          eta_b[which(models=='alpha')]),
                 ddm_params = c(bound,drift),
                 obs=temp_dat,returnFit = F,eta_sep=T,estimate_evidence = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma,targetname = 'fb',fitname = 'cj')
  results_beta_learn <-
    ldc.nn.fit.w(params=c(a0,b0,1,eta_a[which(models=='beta')],
                          eta_b[which(models=='beta')]),
                 ddm_params = c(bound,drift),
                 obs=temp_dat,returnFit = F,eta_sep=T,estimate_evidence = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma,targetname = 'fb',fitname = 'cj')
  results_both_learn <-
    ldc.nn.fit.w(params=c(a0,b0,1,eta_a[which(models=='both')],
                          eta_b[which(models=='both')]),
                 ddm_params = c(bound,drift),
                 obs=temp_dat,returnFit = F,eta_sep=T,estimate_evidence = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma,targetname = 'fb',fitname = 'cj')
  results_no_learn <-
    ldc.nn.fit.w(params=c(a0,b0,1,eta_a[which(models=='no')],
                          eta_b[which(models=='no')]),
                 ddm_params = c(bound,drift),
                 obs=temp_dat,returnFit = F,eta_sep=T,estimate_evidence = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma,targetname = 'fb',fitname = 'cj')
  
  simDat_alpha[simDat_alpha$sub==iter_manip,'cj'] <- results_alpha_learn$pred
  simDat_both[simDat_both$sub==iter_manip,'cj'] <- results_both_learn$pred
  simDat_no[simDat_no$sub==iter_manip,'cj'] <- results_no_learn$pred
  simDat_beta[simDat_beta$sub==iter_manip,'cj'] <- results_beta_learn$pred
}

simDat_no$evidence <- simDat_no$evidence_mean_ev
simDat_alpha$evidence <- simDat_alpha$evidence_mean_ev
simDat_beta$evidence <- simDat_beta$evidence_mean_ev
simDat_both$evidence <- simDat_both$evidence_mean_ev

# write.csv(simDat_alpha,"simDat_alpha_model_recovery_ideal.csv",row.names = F)
# write.csv(simDat_both,"simDat_both_model_recovery_ideal.csv",row.names = F)
# write.csv(simDat_beta,"simDat_beta_model_recovery_ideal.csv",row.names = F)
# write.csv(simDat_no,"simDat_no_model_recovery_ideal.csv",row.names = F)
# Compute rolling mean per subject ----------------------------------------
n <- 25 # Rolling mean window size
n_err <- 25
Ntrials <- ntrial*length(difflevels)
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
                     cor=c(0,1),sub=rep(1:length(fb_manip),each=Ntrials*2))

cj_ma <- merge(trial_conf_sub,trials,all=T)
cj_ma_no <- merge(trial_conf_no,trials,all=T)
cj_ma_alpha <- merge(trial_conf_alpha,trials,all=T)
cj_ma_beta <- merge(trial_conf_beta,trials,all=T)

ma <- function(x,n,names){
  return(rollapply(x[,names], width=n, FUN=function(x) mean(x, na.rm=TRUE),partial=TRUE, align="center"))
}
for (s in 1:length(fb_manip)) {
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

with(simDat_both,aggregate(cj,list(switch,cor,manip),mean))

# Plot data ---------------------------------------------------------------

Nphase_block <- 50
Ntrials_phase <- ntrial_per_switch*length(difflevels)
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

title_line <- -2
cex.title <- 3
cex.lab <- 3
cex.axis <- 2
cex.legend <- 2

go_to('new_sim')
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
     xlab = paste("Consecutive groups of",Ntrials_phase/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
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
     xlab = paste("Consecutive groups of",Ntrials_phase/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
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
     xlab = paste("Consecutive groups of",Ntrials_phase/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
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
     xlab = paste("Consecutive groups of",Ntrials_phase/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
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
     xlab = paste("Consecutive groups of",Ntrials_phase/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
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
     xlab = paste("Consecutive groups of",Ntrials_phase/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
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
     xlab = paste("Consecutive groups of",Ntrials_phase/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
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
     xlab = paste("Consecutive groups of",Ntrials_phase/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
lines(conf_beta$plus_0, type = 'b', pch = 17, col = VERMILLION, lwd = 2, lty = 2)
lines(conf_beta$plus_1, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 1)
lines(conf_beta$minus_1, type = 'b', pch = 16, col = BLUE, lwd = 2, lty = 1)
legend("top",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)
dev.off()
setwd('..')
