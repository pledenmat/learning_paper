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
go_to("data")

Data <- readbulk::read_bulk()

# Let's check the number of rows per participant
table(Data$sub)

# Two participants were given the same number, let's separate them
Data[Data$File=="learn_alfa_sub132_1.csv",'sub'] <- 1320

# Filter participants who finished the experiment
finished_exp <- as.numeric(names(table(Data$sub)[table(Data$sub)==760]))
Data <- subset(Data, sub %in% finished_exp)

questions <- Data[complete.cases(Data$questionID),
                  c("sub","age","gender","handedness","questionID","questionResp")]

Data <- Data[complete.cases(Data$block),]
rm_col <- c("File","questionID","questionResp","X","ï..sub")
Data <- Data[,-which(names(Data) %in% rm_col)]


subs <- unique(Data$sub); Nsub <- length(subs)
Ntrials <- dim(Data)[1]/Nsub

Data$response <- 0
Data[Data$resp=="['n']",'response'] <- 1

## Diagnostic plot per participant and task + chance performance testing
exclusion <- c()
tasks <- unique(Data$task)
par(mfrow=c(2,2))
for(i in 1:Nsub){
  tempDat <- subset(Data,sub==subs[i])
  acc_block <- with(tempDat,aggregate(cor,by=list(block=block),mean))
  bias_block <- with(tempDat,aggregate(response,by=list(block=block),mean))
  if (plots) {
    par(mar=c(3,4,2,0))
    plot(acc_block,ylab="Acc (.) and bias (x)",frame=F,ylim=c(0,1));abline(h=.5,lty=2,col="grey")
    points(bias_block,pch=4)
    plot(tempDat$rt,frame=F,col=c("black"),main=paste('subject',subs[i]),ylab="RT")
    plot(tempDat$cj,frame=F,col=c("black"),ylim=c(1,6),ylab="conf")
    plot(tempDat$RTconf,frame=F,col=c("black"),ylab="RT_conf")
  }
  for (p in unique(Data$phase)) {
    temp <- subset(tempDat,phase==p)
    test <- binom.test(length(temp$cor[temp$cor==1]),n=length(temp$cor),alternative = "greater")
    print(paste("In t0, sub",subs[i],"p =", round(test$p.value,3),"compared to chance"))
    if (test$p.value > .05) {
      exclusion <- c(exclusion,subs[i])
    }
  }
}
par(mar=c(5,4,4,2)+.1,mfrow=c(1,1))

#' Filter out participants who reported only one confidence level 
#' more than 85% of the time
conf_count <- with(Data,aggregate(cj,by=list(sub=sub),max_count))
conf_count$x <- conf_count$x/Ntrials
exclusion <- c(exclusion, unique(conf_count[conf_count$x>.85,"sub"]))

Data <- subset(Data,!(sub %in% exclusion))

age_alpha <- with(Data,aggregate(age,by=list(sub),mean))
summary(age_alpha$x)
gender_alpha <- table(Data$gender)

# Convert RTs to seconds
Data$rt <- Data$rt/1000
Data$RTconf <- Data$RTconf/1000

# Data <- subset(Data,rt>.2 & rt<5) #Trim RTs
# Data <- subset(Data,RTconf<5) #Trim confRTs

Data$response[Data$response==0] <- -1

Data$condition <- as.factor(Data$condition)
Data$difflevel <- as.factor(Data$difflevel)
Data$sub <- as.factor(Data$sub)
subs <- unique(Data$sub); Nsub <- length(subs)

Data[Data$manip=='alfa','manip'] <- 'alpha'
names(Data)[1] <- "trial"

summary(Data)

Data$cj <- Data$cj/6
Data$fb <- Data$fb/100
Data$rt2 <- Data$rt + Data$RTconf

Data$group <- "plus_first"
minus_first <- unique(subset(Data,phase==0&condition=='minus')$sub)
Data[Data$sub %in% minus_first,"group"] <- "minus_first"

# setwd(curdir)
Data_trim <- subset(Data,RTconf<5 & rt<5 & rt>.2)
Data_trim <- subset(Data,RTconf<5 & rt<5 & rt>.2 & trial>4)
# write.csv(Data_trim,file = "alternating_fb_mod_trim_skip.csv",row.names = F)
# 
# Data <- read.csv("alternating_fb_mod.csv")
Data_full <- Data
Data <- Data_trim

Nphase_trial <- length(unique(Data$withinphasetrial))
Nphase_block <- 4
Data$phase_block <-  Data$withinphasetrial %/% (Nphase_trial/Nphase_block)

Data_alpha <- subset(Data,manip=='alpha')
Data_beta <- subset(Data,manip=='beta')
# Behavior analysis -----------------------------------------------------------
if (stat_tests) {
  # RT
  m.int <- lmer(log(rt)~condition*difflevel*manip + (1|sub),data = Data,REML = F)
  m.cond <- lmer(log(rt)~condition*difflevel*manip + (condition|sub),data = Data,REML = F)
  anova(m.int,m.cond)
  m.diff <- lmer(rt~condition*difflevel*manip + (difflevel|sub),data = Data,REML = F,
                 lmerControl(optimizer = 'bobyqa')) #Singular fit
  leveneTest(residuals(m.cond) ~ Data$condition) #Homogeneity of variance
  qqmath(m.cond) #Normality
  anova(m.cond) #Results
  
  # Accuracy
  m.int <- glmer(cor~condition*difflevel + (1|sub),data=Data,family=binomial)
  m.cond <- glmer(cor~condition*difflevel + (condition|sub),data=Data,family=binomial); 
  anova(m.int,m.cond)
  m.diff <- glmer(cor~condition*difflevel*manip + (difflevel|sub),data=Data,family=binomial); #singular
  leveneTest(residuals(m.cond) ~ Data$condition) #Homogeneity of variance
  plot(m.int)
  Anova(m.int)
  
  # Confidence
  Data$cor <- as.factor(Data$cor)
  m.int <- lmer(cj ~ condition*cor*difflevel*manip + (1|sub),data = Data,REML = F); 
  m.cond <- lmer(cj ~ condition*cor*difflevel*manip + (condition|sub),data = Data,REML = F); 
  anova(m.int,m.cond)
  m.cond.acc <- lmer(cj ~ condition*cor*difflevel*manip + (cor + condition|sub),data = Data, REML = F)
  anova(m.cond,m.cond.acc)
  m.all <- lmer(cj ~ condition*cor*difflevel*manip + (cor + condition + difflevel|sub),
                data = Data, REML = F) # Failed to converge
  m.cond.acc.x <- lmer(cj ~ condition*cor*difflevel*manip + (cor * condition|sub),
                       data = Data, REML = F,control = lmerControl(optimizer = 'bobyqa'))
  anova(m.cond.acc,m.cond.acc.x)
  plot(resid(m.cond.acc.x),Data$cj) #Linearity
  plot(m.cond.acc.x)
  leveneTest(residuals(m.cond.acc.x) ~ Data$cor*Data$condition*Data$difflevel*Data$manip) #Homogeneity of variance
  qqmath(m.cond.acc.x) #Normality
  anova(m.cond.acc.x) #Results
  
  # No condition*manip interaction but let's still check each parameter separately
  # Alpha
  m.int <- lmer(cj ~ condition*cor*difflevel + (1|sub),data = Data_alpha,REML = F); 
  m.cond <- lmer(cj ~ condition*cor*difflevel + (condition|sub),data = Data_alpha,REML = F); 
  anova(m.int,m.cond)
  m.cond.acc <- lmer(cj ~ condition*cor*difflevel + (cor + condition|sub),data = Data_alpha, REML = F)
  anova(m.cond,m.cond.acc)
  m.cond.acc.x <- lmer(cj ~ condition*cor*difflevel + (cor * condition|sub),
                       data = Data_alpha, REML = F,control = lmerControl(optimizer = 'bobyqa'))
  anova(m.cond.acc,m.cond.acc.x)
  plot(resid(m.cond.acc.x),Data_alpha$cj) #Linearity
  leveneTest(residuals(m.cond.acc.x) ~ Data_alpha$cor*Data_alpha$condition*Data_alpha$difflevel) #Homogeneity of variance
  qqmath(m.cond.acc.x) #Normality
  anova(m.cond.acc.x) #Results
  
  # Beta
  m.int <- lmer(cj ~ condition*cor*difflevel + (1|sub),data = Data_beta,REML = F); 
  m.cond <- lmer(cj ~ condition*cor*difflevel + (condition|sub),data = Data_beta,REML = F); 
  anova(m.int,m.cond)
  m.cond.acc <- lmer(cj ~ condition*cor*difflevel + (cor + condition|sub),data = Data_beta, REML = F)
  anova(m.cond,m.cond.acc)
  m.cond.acc.x <- lmer(cj ~ condition*cor*difflevel + (cor * condition|sub),
                       data = Data_beta, REML = F,control = lmerControl(optimizer = 'bobyqa'))
  anova(m.cond.acc,m.cond.acc.x)
  plot(resid(m.cond.acc.x),Data_beta$cj) #Linearity
  leveneTest(residuals(m.cond.acc.x) ~ Data_beta$cor*Data_beta$condition*Data_beta$difflevel) #Homogeneity of variance
  qqmath(m.cond.acc.x) #Normality
  anova(m.cond.acc.x) #Results
  
  Data$cor <- as.numeric(Data$cor)
}
# Retrieve fits -----------------------------------------------------------
beta_input <- .1
Nupdate_per_trial <- 1
model <- "allpar"
dt <- .001; sigma <- .1
binning <- F

conditions <- sort(unique(Data$condition))
difflevels <- sort(unique(Data$difflevel)) # Important to sort to match with drift order

ntrial <- Ntrials
par_trace <- data.frame(trial=0:(ntrial-1),sub=rep(subs,each=ntrial*2),
                        condition=rep(c("plus","minus"),each=ntrial),
                        alpha=NA,beta=NA)


totlen <- length(conditions)*length(difflevels)*length(subs)
par <- data.frame(bound = NA, drift = NA, ter = NA, eta = NA,
                  cost_ddm = NA, cost_ldc = NA,
                  sub = rep(subs, each = length(conditions)*length(difflevels)),
                  condition = rep(conditions,each = length(difflevels),
                                  length.out=totlen),
                  difflevel = rep(difflevels, length.out = totlen),
                  manip=NA)
go_to("fit")
go_to("alternating_fb")
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
  
  par[par$sub==subs[s],"bound"] <- ddm.results$optim$bestmem[1]
  par[par$sub==subs[s],"ter"] <- ddm.results$optim$bestmem[2]
  par[par$sub==subs[s],"z"] <- ddm.results$optim$bestmem[3]
  par[par$sub==subs[s],"vratio"] <- ddm.results$optim$bestmem[4]
  par[par$sub==subs[s]&par$difflevel==difflevels[1],"drift"] <- ddm.results$optim$bestmem[5]
  par[par$sub==subs[s]&par$difflevel==difflevels[2],"drift"] <- ddm.results$optim$bestmem[6]
  par[par$sub==subs[s]&par$difflevel==difflevels[3],"drift"] <- ddm.results$optim$bestmem[7]
  par[par$sub==subs[s],"cost_ddm"] <- ddm.results$optim$bestval
  

  ldc_file <- paste0('ldc_nn/trim_skip/batch_',Nupdate_per_trial,'/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
  }else{
    ldc_save <- paste0('ldc_nn/trim_skip/batch_',Nupdate_per_trial,'/slurm-55185173_',s,'.out')
    test <- readLines(ldc_save)
    check_finish <- "***** summary of DEoptim object ***** "
    if (test[length(test)-5] == check_finish) {
      fit_par <- unlist(strsplit(test[length(test)-4],split=" "))
      a0 <- as.numeric(fit_par[7])
      b0 <- as.numeric(fit_par[8])
      eta <- as.numeric(fit_par[10])
      cost <- as.numeric(unlist(strsplit(test[length(test)-3],split=" "))[8])
      ldc.results <- list(optim=list(bestmem=c(a0,b0,1,eta),bestval=cost))
      save(ldc.results,file=ldc_file)
    }
  }

  par[par$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
  par[par$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
  par[par$sub==subs[s],"eta"] <- ldc.results$optim$bestmem[4]
  par[par$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  
  results <-
    ldc.nn.fit.w(params=c(mean(par[par$sub==subs[s],"a0"]),
                          mean(par[par$sub==subs[s],"b0"]),1,
                          mean(par[par$sub==subs[s],"eta"])),
                 ddm_params = ddm_params,
                 obs=temp_dat,returnFit = F,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma)
  Data[Data$sub==subs[s],'cj_pred'] <- results$pred

  par_trace[par_trace$sub==subs[s],"alpha"] <-
    c(results$trace[,1],
      rep(NA,ntrial-length(results$trace[seq(Nupdate_per_trial,nrow(results$trace),Nupdate_per_trial),1])))
  par_trace[par_trace$sub==subs[s],"beta"] <-
    c(results$trace[,2],
      rep(NA,ntrial-length(results$trace[seq(Nupdate_per_trial,nrow(results$trace),Nupdate_per_trial),2])))


  par[par$sub==subs[s],"manip"] <- unique(temp_dat$manip)
}
summary(par)
par <- par[complete.cases(par$bound),]

par_trace <- par_trace[complete.cases(par_trace),]
par_trace <- subset(par_trace,condition=='plus')
par_trace <- par_trace[order(par_trace$sub),]

Data$alpha <- par_trace$alpha
Data$beta <- par_trace$beta
Data_alpha <- subset(Data,manip=='alpha')
Data_beta <- subset(Data,manip=='beta')

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
go_to("trim_skip")

jpeg(filename = "traces.jpg",units = 'cm',width = 42,height = 30,res=300)
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
# RTconf_rm <- subset(Data,trial>5&RTconf>5)
# RTconf_rm <- unique(RTconf_rm$sub)
# Data_beta <- subset(Data_beta,!(sub %in% RTconf_rm))
# Data_alpha <- subset(Data_alpha,!(sub %in% RTconf_rm))

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





# Trace per sub --------------------------------------------------------
go_to("trace_per_sub")

for (s in subs) {
  temp_ma <- subset(cj_ma,sub==s&cor==1)
  temp_pred_ma <- subset(cj_pred_ma,sub==s&cor==1)
  switches <- which(subset(Data,sub==s)$trial %in% 
                      seq(Ntrials_phase,Ntrials-1,Ntrials_phase))
  jpeg(paste0("trace_",s,".jpg"),width=width,height=height,units = 'cm',res=300)
  plot(temp_ma$cj,type='l',col=VERMILLION,ylim=c(.5,1),
       main= paste("Confidence in correct trials, sub",s),
       xlab = "Trial", ylab = "Confidence")
  lines(temp_pred_ma$cj,col=VERMILLION,lty=2)
  abline(v=switches,lty=2,col='lightgrey')
  legend("bottomleft",horiz=T,lty=c(1,2),legend=c("Behaviour","Model"),bty='n')
  dev.off()
}
setwd("..")
# Testing wild mixed model ------------------------------------------------

m.int <- lmer(data = Data, cj ~ condition*cor*difflevel*group*manip + (1|sub), REML = F)
m.cond <- lmer(data = Data, cj ~ condition*cor*difflevel*group*manip + (condition|sub), REML = F)
anova(m.int,m.cond)
m.cond.cor <- lmer(data = Data, cj ~ condition*cor*difflevel*group*manip + (condition + cor|sub), REML = F)
anova(m.cond,m.cond.cor)
# m.cond.cor.diff <- lmer(data = Data, cj ~ condition*cor*difflevel*group*manip + 
#                           (condition + cor + difflevel|sub), REML = F) # Singular
anova(m.cond.cor)
emm <- emmeans(m.cond.cor, ~ condition | group)
pairs(emm)
emm <- emmeans(m.cond.cor, ~ condition | manip)
pairs(emm)
emm <- emmeans(m.cond.cor, ~ condition | cor | manip*group)
pairs(emm)

conf_group <- with(Data,aggregate(cj,by=list(group,manip,condition),mean))
conf_group_sd <- with(Data,aggregate(cj,by=list(group,manip,condition),sd))
conf_group_length <- with(Data,aggregate(cj,by=list(group,manip,condition),length))
names(conf_group) <- c('group','manip','condition','cj')
names(conf_group_sd) <- c('group','manip','condition','cj')
names(conf_group_length) <- c('group','manip','condition','Nsub')
conf_group_length$Nsub <- conf_group_length$Nsub/Ntrials*2
conf_group_sd$cj <- conf_group_sd$cj/sqrt(conf_group_length$Nsub)

conf_alpha <- as.matrix(cast(subset(conf_group,manip=='alpha'),condition~group))
conf_beta <- as.matrix(cast(subset(conf_group,manip=='beta'),condition~group))
par(mfrow=c(1,2))
barplot(conf_alpha, main='Alpha experiment',ylim=c(0,1),
        col=colors()[c(23,89)] , 
        border="white", 
        font.axis=2, 
        beside=T, 
        legend=rownames(conf_alpha),
        args.legend = list(x='top',bty='n',horiz=T),
        xlab="group", 
        font.lab=2)
barplot(conf_beta, main='Beta experiment',ylim=c(0,1),
        col=colors()[c(23,89)] , 
        border="white", 
        font.axis=2, 
        beside=T, 
        # legend=rownames(conf_beta),
        xlab="group", 
        font.lab=2)
par(mfrow=c(1,1))

# Feedback presented ExpA ------------------------------------------------
N_temp <- length(unique(Data_alpha$sub))

fbminus <- with(subset(Data_alpha,condition=="minus"),aggregate(fb,by=list(sub,difflevel,cor),mean));
names(fbminus) <- c('sub','difflevel','cor','fb')
fbminus_cor <- subset(fbminus,cor==1); fbminus_err <- subset(fbminus,cor==0)
fbminus_cor <- cast(fbminus_cor,sub~difflevel); fbminus_err <- cast(fbminus_err,sub~difflevel)

fbplus <- with(subset(Data_alpha,condition=="plus"),aggregate(fb,by=list(sub,difflevel,cor),mean));
names(fbplus) <- c('sub','difflevel','cor','fb')
fbplus_cor <- subset(fbplus,cor==1); fbplus_err <- subset(fbplus,cor==0)
fbplus_cor <- cast(fbplus_cor,sub~difflevel); fbplus_err <- cast(fbplus_err,sub~difflevel)


# Drop subject column
xminus_cor <- fbminus_cor[,c(2:4)];xplus_cor <- fbplus_cor[,c(2:4)]
xminus_err <- fbminus_err[,c(2:4)];xplus_err <- fbplus_err[,c(2:4)]
n <- length(xminus_err)

xminus_cor <- xminus_cor[,c("hard","medium","easy")];
xplus_cor <- xplus_cor[,c("hard","medium","easy")]
xminus_err <- xminus_err[,c("hard","medium","easy")];
xplus_err <- xplus_err[,c("hard","medium","easy")]

stripchart(xminus_cor,ylim=c(0,100), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           yaxt = 'n',xlab="",ylab = "")
title(ylab = "Feedback", xlab = "Trial difficulty", line = 2.5,cex.lab=cex.lab.rel)
axis(1,at=0:(n-1),labels=names(xminus_cor), cex.axis=cex.axis);
axis(2, seq(0,100,20), cex.axis=cex.axis)
means <- sapply(xminus_cor, mean)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_minus,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xminus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_minus)
means <- sapply(xplus_cor, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_plus,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xplus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_plus)

means <- sapply(xminus_err, mean, na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_minus,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xminus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_minus)
means <- sapply(xplus_err, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_plus,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xplus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_plus)

legend("bottomleft",legend=c("Correct trials","Error trials"),
       title = NULL,lty=c("dashed","dotted"),bty = "n",inset=0,
       cex = cex.legend, lwd=lwd.dat,seg.len=1.5)
legend("bottomright",legend=c("Minus","Plus"),
       title = NULL,pch=rep(16,3),bty = "n",inset=0,
       cex = cex.legend,col=c(col_minus,col_plus))
# Feedback presented ExpB ------------------------------------------------
N_temp <- length(unique(Data_beta$sub))

fbminus <- with(subset(Data_beta,condition=="minus"),aggregate(fb,by=list(sub,difflevel,cor),mean));
names(fbminus) <- c('sub','difflevel','cor','fb')
fbminus_cor <- subset(fbminus,cor==1); fbminus_err <- subset(fbminus,cor==0)
fbminus_cor <- cast(fbminus_cor,sub~difflevel); fbminus_err <- cast(fbminus_err,sub~difflevel)

fbplus <- with(subset(Data_beta,condition=="plus"),aggregate(fb,by=list(sub,difflevel,cor),mean));
names(fbplus) <- c('sub','difflevel','cor','fb')
fbplus_cor <- subset(fbplus,cor==1); fbplus_err <- subset(fbplus,cor==0)
fbplus_cor <- cast(fbplus_cor,sub~difflevel); fbplus_err <- cast(fbplus_err,sub~difflevel)


# Drop subject column
xminus_cor <- fbminus_cor[,c(2:4)];xplus_cor <- fbplus_cor[,c(2:4)]
xminus_err <- fbminus_err[,c(2:4)];xplus_err <- fbplus_err[,c(2:4)]
n <- length(xminus_err)

xminus_cor <- xminus_cor[,c("hard","medium","easy")];
xplus_cor <- xplus_cor[,c("hard","medium","easy")]
xminus_err <- xminus_err[,c("hard","medium","easy")];
xplus_err <- xplus_err[,c("hard","medium","easy")]

stripchart(xminus_cor,ylim=c(0,100), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           yaxt = 'n',xlab="",ylab = "")
title(ylab = "Feedback", xlab = "Trial difficulty", line = 2.5,cex.lab=cex.lab.rel)
axis(1,at=0:(n-1),labels=names(xminus_cor), cex.axis=cex.axis);
axis(2, seq(0,100,20), cex.axis=cex.axis)
means <- sapply(xminus_cor, mean)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_minus,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xminus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_minus)
means <- sapply(xplus_cor, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_plus,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xplus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_plus)

means <- sapply(xminus_err, mean, na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_minus,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xminus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_minus)
means <- sapply(xplus_err, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_plus,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xplus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_plus)

legend("bottomleft",legend=c("Correct trials","Error trials"),
       title = NULL,lty=c("dashed","dotted"),bty = "n",inset=0,
       cex = cex.legend, lwd=lwd.dat,seg.len=1.5)



# Plot correlation between observed and predicted -------------------------
plot(Data$cj,Data$cj_pred)
library(vioplot)
cj_lab <- sort(unique(Data$cj))
jpeg(filename="scatter_pred_behav.jpg",width=16,height=8,units='cm',res=300)
with(Data , vioplot( 
  cj_pred[cj==cj_lab[1]] , cj_pred[cj==cj_lab[2]], cj_pred[cj==cj_lab[3]],  
  cj_pred[cj==cj_lab[4]] , cj_pred[cj==cj_lab[5]], cj_pred[cj==cj_lab[6]],  
  col=rgb(0.1,0.4,0.7,0.7) , names=c(1,2,3,4,5,6),
  xlab = "Confidence (behavior)", ylab = "Confidence (model)",
  main = paste("cor =",round(cor(Data$cj,Data$cj_pred),3))
))

jpeg(filename="density_pred_behav.jpg",width=16,height=12,units='cm',res=300)
x_coord <- 0:5
plot.new()
plot.window(xlim=c(-.1,6), ylim=c(-.1,1),xaxt='n',
            xlab = "Confidence (behavior)", ylab = "Confidence (model)",
            main=paste("cor =",round(cor(Data$cj,Data$cj_pred),3)))
# plot(a_free$`0`,, col="white",frame=F,, 
#      ylab="",yaxt='n',xlab="");
axis(2,at=seq(0,1,.2),label=seq(0,1,.2))
axis(1,at=x_coord,label=cj_lab*6)
mtext(text = "Confidence (model)",side = 2, line = 2.5)
mtext(text = "Confidence (behavior)",side = 1, line = 2.5)
mtext(text = paste("cor =",round(cor(Data$cj,Data$cj_pred),3)),side = 3,line=1,
      cex = 1.5)
# segments(x0=0,y0=0,y1=2,lty=2,col="lightgrey")
for(i in 1:length(x_coord)){
  polygon((rescale(density(subset(Data,cj==cj_lab[i])$cj_pred)$y,to=c(0,.9)))+x_coord[i],
                  density(subset(Data,cj==cj_lab[i])$cj_pred)$x,
          col=rgb(0.1,0.4,0.7,0.7),border=F)
}
abline(lm(Data$cj~Data$cj_pred)$coef[1],lm(Data$cj~Data$cj_pred)$coef[2]/6)
dev.off()

# Split trials in each phase in 4 -----------------------------------------
m.int <- lmer(data = Data_alpha, cj ~ condition*phase_block*cor + (1|sub),REML = F)
m.cond <- lmer(data = Data_alpha, cj ~ condition*phase_block*cor + (condition|sub),REML = F)
anova(m.int,m.cond)
m.cond.cor <- lmer(data = Data_alpha, cj ~ condition*phase_block*cor + (condition+cor|sub),REML = F)
anova(m.cond,m.cond.cor)
anova(m.cond.cor)
emm <- emmeans(m.cond.cor, ~ condition | phase_block)
pairs(emm)

m.int <- lmer(data = Data_beta, cj ~ condition*phase_block*cor + (1|sub),REML = F)
m.cond <- lmer(data = Data_beta, cj ~ condition*phase_block*cor + (condition|sub),REML = F)
anova(m.int,m.cond)
m.cond.cor <- lmer(data = Data_beta, cj ~ condition*phase_block*cor + (condition+cor|sub),REML = F)
anova(m.cond,m.cond.cor)
anova(m.cond.cor)
emm <- emmeans(m.cond.cor, ~ condition | cor | phase_block)
pairs(emm)

Nalpha <- length(unique(Data_alpha$sub))
Nbeta <- length(unique(Data_beta$sub))
conf_group <- with(Data,aggregate(cj,by=list(phase_block,manip,condition),mean))
conf_group_sd <- with(Data,aggregate(cj,by=list(phase_block,manip,condition),sd))
names(conf_group) <- c('phase_block','manip','condition','cj')
names(conf_group_sd) <- c('phase_block','manip','condition','cj')

conf_alpha <- cast(subset(conf_group,manip=='alpha'),phase_block~condition)
conf_beta <- cast(subset(conf_group,manip=='beta'),phase_block~condition)
conf_alpha_sd <- subset(conf_group_sd,manip=='alpha')
conf_alpha_sd$cj <- conf_alpha_sd$cj/sqrt(Nalpha)
conf_alpha_sd <- cast(conf_alpha_sd, phase_block ~ condition)
conf_beta_sd <- subset(conf_group_sd,manip=='beta')
conf_beta_sd$cj <- conf_beta_sd$cj/sqrt(Nbeta)
conf_beta_sd <- cast(conf_beta_sd, phase_block ~ condition)

jpeg('trace_aggreg.jpg', width = 20, height = 12, units = 'cm', res = 300)
par(mfrow=c(1,2))
plot(conf_alpha$minus,ylim=c(.79,.88),main='Alpha experiment',col = BLUE, type = 'b',
     lty = 1, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",
     xlab = "Within phase block of trials")
axis(1, at = 1:4, labels = 1:4)
lines(conf_alpha$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2)
error.bar(1:length(conf_alpha$minus),conf_alpha$minus,conf_alpha_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:length(conf_alpha$plus),conf_alpha$plus,conf_alpha_sd$plus,
          lwd=2, col = VERMILLION)
legend("bottom",legend = c('Plus','Minus'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n')


plot(conf_beta$minus,ylim=c(.79,.88),main='Beta experiment',col = BLUE, type = 'b',
     lty = 1, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",
     xlab = "Within phase block of trials")
axis(1, at = 1:4, labels = 1:4)
lines(conf_beta$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2)
error.bar(1:length(conf_beta$minus),conf_beta$minus,conf_beta_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:length(conf_beta$plus),conf_beta$plus,conf_beta_sd$plus,
          lwd=2, col = VERMILLION)
par(mfrow=c(1,1))
dev.off()

# Aggregating behavior
Nalpha <- length(unique(Data_alpha$sub))
Nbeta <- length(unique(Data_beta$sub))
conf_group <- with(Data,aggregate(cj,by=list(phase_block,manip,condition,cor),mean))
conf_group_sd <- with(Data,aggregate(cj,by=list(phase_block,manip,condition,cor),sd))
names(conf_group) <- c('phase_block','manip','condition','cor','cj')
names(conf_group_sd) <- c('phase_block','manip','condition','cor','cj')

conf_alpha <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor)
conf_beta <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor)
conf_alpha_sd <- subset(conf_group_sd,manip=='alpha')
conf_alpha_sd$cj <- conf_alpha_sd$cj/sqrt(Nalpha)
conf_alpha_sd <- cast(conf_alpha_sd, phase_block ~ condition + cor)
conf_beta_sd <- subset(conf_group_sd,manip=='beta')
conf_beta_sd$cj <- conf_beta_sd$cj/sqrt(Nbeta)
conf_beta_sd <- cast(conf_beta_sd, phase_block ~ condition + cor)

# Aggregating model fits
conf_group_pred <- with(Data,aggregate(cj_pred,by=list(phase_block,manip,condition,cor),mean))
conf_group_pred_sd <- with(Data,aggregate(cj_pred,by=list(phase_block,manip,condition,cor),sd))
names(conf_group_pred) <- c('phase_block','manip','condition','cor','cj')
names(conf_group_pred_sd) <- c('phase_block','manip','condition','cor','cj')

conf_alpha_pred <- cast(subset(conf_group_pred,manip=='alpha'),phase_block~condition+cor)
conf_beta_pred <- cast(subset(conf_group_pred,manip=='beta'),phase_block~condition+cor)
conf_alpha_pred_sd <- subset(conf_group_pred_sd,manip=='alpha')
conf_alpha_pred_sd$cj <- conf_alpha_pred_sd$cj/sqrt(Nalpha)
conf_alpha_pred_sd <- cast(conf_alpha_pred_sd, phase_block ~ condition + cor)
conf_beta_pred_sd <- subset(conf_group_pred_sd,manip=='beta')
conf_beta_pred_sd$cj <- conf_beta_pred_sd$cj/sqrt(Nbeta)
conf_beta_pred_sd <- cast(conf_beta_pred_sd, phase_block ~ condition + cor)

xlen <- nrow(conf_alpha)
jpeg('trace_aggreg_cor.jpg', width = 20, height = 12, units = 'cm', res = 300)
par(mfrow=c(1,2))
plot(conf_alpha$minus_0,ylim=c(.7,.9),main='Alpha experiment',col = BLUE, type = 'b',
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",
     xlab = "Within phase block of trials")
axis(1, at = 1:4, labels = 1:4)
lines(conf_alpha$plus_0, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
lines(conf_alpha$plus_1, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 1)
lines(conf_alpha$minus_1, type = 'b', pch = 16, col = BLUE, lwd = 2, lty = 1)
error.bar(1:xlen,conf_alpha$minus_0,conf_alpha_sd$minus_0,
          lwd=2, col = BLUE)
error.bar(1:xlen,conf_alpha$plus_0,conf_alpha_sd$plus_0,
          lwd=2, col = VERMILLION)
error.bar(1:xlen,conf_alpha$minus_1,conf_alpha_sd$minus_1,
          lwd=2, col = BLUE)
error.bar(1:xlen,conf_alpha$plus_1,conf_alpha_sd$plus_1,
          lwd=2, col = VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_alpha_pred$minus_0 + conf_alpha_pred_sd$minus_0,
        (conf_alpha_pred$minus_0 - conf_alpha_pred_sd$minus_0)[xlen:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:xlen,xlen:1),c(conf_alpha_pred$minus_1 + conf_alpha_pred_sd$minus_1,
                           (conf_alpha_pred$minus_1 - conf_alpha_pred_sd$minus_1)[xlen:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:xlen,xlen:1),c(conf_alpha_pred$plus_0 + conf_alpha_pred_sd$plus_0,
                           (conf_alpha_pred$plus_0 - conf_alpha_pred_sd$plus_0)[xlen:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))
polygon(c(1:xlen,xlen:1),c(conf_alpha_pred$plus_1 + conf_alpha_pred_sd$plus_1,
                           (conf_alpha_pred$plus_1 - conf_alpha_pred_sd$plus_1)[xlen:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))
 legend("bottom",legend = c('Plus','Minus'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n')


plot(conf_beta$minus_0,ylim=c(.7,.9),main='Beta experiment',col = BLUE, type = 'b',
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",
     xlab = "Within phase block of trials")
axis(1, at = 1:4, labels = 1:4)
lines(conf_beta$plus_0, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
lines(conf_beta$plus_1, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 1)
lines(conf_beta$minus_1, type = 'b', pch = 16, col = BLUE, lwd = 2, lty = 1)
error.bar(1:length(conf_beta$minus_0),conf_beta$minus_0,conf_beta_sd$minus_0,
          lwd=2, col = BLUE)
error.bar(1:length(conf_beta$plus_0),conf_beta$plus_0,conf_beta_sd$plus_0,
          lwd=2, col = VERMILLION)
error.bar(1:length(conf_beta$minus_1),conf_beta$minus_1,conf_beta_sd$minus_1,
          lwd=2, col = BLUE)
error.bar(1:length(conf_beta$plus_1),conf_beta$plus_1,conf_beta_sd$plus_1,
          lwd=2, col = VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_beta_pred$minus_0 + conf_beta_pred_sd$minus_0,
                           (conf_beta_pred$minus_0 - conf_beta_pred_sd$minus_0)[xlen:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:xlen,xlen:1),c(conf_beta_pred$minus_1 + conf_beta_pred_sd$minus_1,
                           (conf_beta_pred$minus_1 - conf_beta_pred_sd$minus_1)[xlen:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:xlen,xlen:1),c(conf_beta_pred$plus_0 + conf_beta_pred_sd$plus_0,
                           (conf_beta_pred$plus_0 - conf_beta_pred_sd$plus_0)[xlen:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))
polygon(c(1:xlen,xlen:1),c(conf_beta_pred$plus_1 + conf_beta_pred_sd$plus_1,
                           (conf_beta_pred$plus_1 - conf_beta_pred_sd$plus_1)[xlen:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))
par(mfrow=c(1,1))
dev.off()

# Parameter trace aggregated over switches ----------------------------------------------
# Analyze alpha experiment
m.alpha.a <- lmer(data = Data_alpha, alpha ~ phase_block*condition + (phase_block+condition|sub), REML = F,
           control = lmerControl(optimizer = "bobyqa"))
anova(m.alpha.a)
emm <- emmeans(m.alpha.a, ~ condition | phase_block)
pairs(emm)
emm <- emmeans(m.alpha.a, ~ phase_block | condition)
pairs(emm)

m.alpha.b <- lmer(data = Data_alpha, beta ~ phase_block*condition + (condition + phase_block|sub), REML = F,
           control = lmerControl(optimizer = "bobyqa"))
anova(m.alpha.b)
emm <- emmeans(m.alpha.b, ~ phase_block | condition)
pairs(emm)
emm <- emmeans(m.alpha.b, ~ condition | phase_block)
pairs(emm)

# Do the same with beta experiment
m.beta.a <- lmer(data = Data_beta, alpha ~ phase_block*condition + (phase_block+condition|sub), REML = F,
                 control = lmerControl(optimizer = "bobyqa"))
anova(m.beta.a)
emm <- emmeans(m.beta.a, ~ phase_block)
pairs(emm)
emm <- emmeans(m.beta.a, ~ condition | phase_block)
pairs(emm)
emm <- emmeans(m.beta.a, ~ phase_block | condition)
pairs(emm)

m.beta.b <- lmer(data = Data_beta, beta ~ phase_block*condition + (condition + phase_block|sub), REML = F,
                 control = lmerControl(optimizer = "bobyqa"))
anova(m.beta.b)

# Plot all this
Nalpha <- length(unique(Data_alpha$sub))
Nbeta <- length(unique(Data_beta$sub))
conf_group <- with(Data,aggregate(alpha,by=list(phase_block,manip,condition),mean))
conf_group_sd <- with(Data,aggregate(alpha,by=list(phase_block,manip,condition),sd))
names(conf_group) <- c('phase_block','manip','condition','cj')
names(conf_group_sd) <- c('phase_block','manip','condition','cj')

conf_alpha <- cast(subset(conf_group,manip=='alpha'),phase_block~condition)
conf_beta <- cast(subset(conf_group,manip=='beta'),phase_block~condition)
conf_alpha_sd <- subset(conf_group_sd,manip=='alpha')
conf_alpha_sd$cj <- conf_alpha_sd$cj/sqrt(Nalpha)
conf_alpha_sd <- cast(conf_alpha_sd, phase_block ~ condition)
conf_beta_sd <- subset(conf_group_sd,manip=='beta')
conf_beta_sd$cj <- conf_beta_sd$cj/sqrt(Nbeta)
conf_beta_sd <- cast(conf_beta_sd, phase_block ~ condition)

jpeg('alpha_trace_aggreg.jpg', width = 20, height = 12, units = 'cm', res = 300)
par(mfrow=c(1,2))
plot(conf_alpha$minus,ylim=c(0,15),main='Alpha experiment',col = BLUE, type = 'b',
     lty = 1, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Alpha",
     xlab = "Within phase block of trials")
axis(1, at = 1:4, labels = 1:4)
lines(conf_alpha$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2)
error.bar(1:length(conf_alpha$minus),conf_alpha$minus,conf_alpha_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:length(conf_alpha$plus),conf_alpha$plus,conf_alpha_sd$plus,
          lwd=2, col = VERMILLION)
legend("bottom",legend = c('Plus','Minus'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n')


plot(conf_beta$minus,ylim=c(0,15),main='Beta experiment',col = BLUE, type = 'b',
     lty = 1, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Alpha",
     xlab = "Within phase block of trials")
axis(1, at = 1:4, labels = 1:4)
lines(conf_beta$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2)
error.bar(1:length(conf_beta$minus),conf_beta$minus,conf_beta_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:length(conf_beta$plus),conf_beta$plus,conf_beta_sd$plus,
          lwd=2, col = VERMILLION)
par(mfrow=c(1,1))
dev.off()

Nalpha <- length(unique(Data_alpha$sub))
Nbeta <- length(unique(Data_beta$sub))
conf_group <- with(Data,aggregate(beta,by=list(phase_block,manip,condition),mean))
conf_group_sd <- with(Data,aggregate(beta,by=list(phase_block,manip,condition),sd))
names(conf_group) <- c('phase_block','manip','condition','cj')
names(conf_group_sd) <- c('phase_block','manip','condition','cj')

conf_alpha <- cast(subset(conf_group,manip=='alpha'),phase_block~condition)
conf_beta <- cast(subset(conf_group,manip=='beta'),phase_block~condition)
conf_alpha_sd <- subset(conf_group_sd,manip=='alpha')
conf_alpha_sd$cj <- conf_alpha_sd$cj/sqrt(Nalpha)
conf_alpha_sd <- cast(conf_alpha_sd, phase_block ~ condition)
conf_beta_sd <- subset(conf_group_sd,manip=='beta')
conf_beta_sd$cj <- conf_beta_sd$cj/sqrt(Nbeta)
conf_beta_sd <- cast(conf_beta_sd, phase_block ~ condition)

jpeg('beta_trace_aggreg.jpg', width = 20, height = 12, units = 'cm', res = 300)
par(mfrow=c(1,2))
plot(conf_alpha$minus,ylim=c(15,22),main='Alpha experiment',col = BLUE, type = 'b',
     lty = 1, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Beta",
     xlab = "Within phase block of trials")
axis(1, at = 1:4, labels = 1:4)
lines(conf_alpha$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2)
error.bar(1:length(conf_alpha$minus),conf_alpha$minus,conf_alpha_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:length(conf_alpha$plus),conf_alpha$plus,conf_alpha_sd$plus,
          lwd=2, col = VERMILLION)
legend("bottom",legend = c('Plus','Minus'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n')


plot(conf_beta$minus,ylim=c(15,22),main='Beta experiment',col = BLUE, type = 'b',
     lty = 1, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Beta",
     xlab = "Within phase block of trials")
axis(1, at = 1:4, labels = 1:4)
lines(conf_beta$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2)
error.bar(1:length(conf_beta$minus),conf_beta$minus,conf_beta_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:length(conf_beta$plus),conf_beta$plus,conf_beta_sd$plus,
          lwd=2, col = VERMILLION)
par(mfrow=c(1,1))
dev.off()


# Plot parameter aggreg trace with P1 skipped -----------------------------

Data_save <- Data
Data <- subset(Data,phase>0)
Nalpha <- length(unique(Data_alpha$sub))
Nbeta <- length(unique(Data_beta$sub))
conf_group <- with(Data,aggregate(alpha,by=list(phase_block,manip,condition),mean))
conf_group_sd <- with(Data,aggregate(alpha,by=list(phase_block,manip,condition),sd))
names(conf_group) <- c('phase_block','manip','condition','cj')
names(conf_group_sd) <- c('phase_block','manip','condition','cj')

conf_alpha <- cast(subset(conf_group,manip=='alpha'),phase_block~condition)
conf_beta <- cast(subset(conf_group,manip=='beta'),phase_block~condition)
conf_alpha_sd <- subset(conf_group_sd,manip=='alpha')
conf_alpha_sd$cj <- conf_alpha_sd$cj/sqrt(Nalpha)
conf_alpha_sd <- cast(conf_alpha_sd, phase_block ~ condition)
conf_beta_sd <- subset(conf_group_sd,manip=='beta')
conf_beta_sd$cj <- conf_beta_sd$cj/sqrt(Nbeta)
conf_beta_sd <- cast(conf_beta_sd, phase_block ~ condition)

jpeg('alpha_trace_aggreg_p1skip.jpg', width = 20, height = 12, units = 'cm', res = 300)
par(mfrow=c(1,2))
plot(conf_alpha$minus,ylim=c(0,15),main='Alpha experiment',col = BLUE, type = 'b',
     lty = 1, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Alpha",
     xlab = "Within phase block of trials")
axis(1, at = 1:4, labels = 1:4)
lines(conf_alpha$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2)
error.bar(1:length(conf_alpha$minus),conf_alpha$minus,conf_alpha_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:length(conf_alpha$plus),conf_alpha$plus,conf_alpha_sd$plus,
          lwd=2, col = VERMILLION)
legend("bottom",legend = c('Plus','Minus'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n')


plot(conf_beta$minus,ylim=c(0,15),main='Beta experiment',col = BLUE, type = 'b',
     lty = 1, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Alpha",
     xlab = "Within phase block of trials")
axis(1, at = 1:4, labels = 1:4)
lines(conf_beta$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2)
error.bar(1:length(conf_beta$minus),conf_beta$minus,conf_beta_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:length(conf_beta$plus),conf_beta$plus,conf_beta_sd$plus,
          lwd=2, col = VERMILLION)
par(mfrow=c(1,1))
dev.off()

Nalpha <- length(unique(Data_alpha$sub))
Nbeta <- length(unique(Data_beta$sub))
conf_group <- with(Data,aggregate(beta,by=list(phase_block,manip,condition),mean))
conf_group_sd <- with(Data,aggregate(beta,by=list(phase_block,manip,condition),sd))
names(conf_group) <- c('phase_block','manip','condition','cj')
names(conf_group_sd) <- c('phase_block','manip','condition','cj')

conf_alpha <- cast(subset(conf_group,manip=='alpha'),phase_block~condition)
conf_beta <- cast(subset(conf_group,manip=='beta'),phase_block~condition)
conf_alpha_sd <- subset(conf_group_sd,manip=='alpha')
conf_alpha_sd$cj <- conf_alpha_sd$cj/sqrt(Nalpha)
conf_alpha_sd <- cast(conf_alpha_sd, phase_block ~ condition)
conf_beta_sd <- subset(conf_group_sd,manip=='beta')
conf_beta_sd$cj <- conf_beta_sd$cj/sqrt(Nbeta)
conf_beta_sd <- cast(conf_beta_sd, phase_block ~ condition)

jpeg('beta_trace_aggreg_p1skip.jpg', width = 20, height = 12, units = 'cm', res = 300)
par(mfrow=c(1,2))
plot(conf_alpha$minus,ylim=c(15,22),main='Alpha experiment',col = BLUE, type = 'b',
     lty = 1, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Beta",
     xlab = "Within phase block of trials")
axis(1, at = 1:4, labels = 1:4)
lines(conf_alpha$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2)
error.bar(1:length(conf_alpha$minus),conf_alpha$minus,conf_alpha_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:length(conf_alpha$plus),conf_alpha$plus,conf_alpha_sd$plus,
          lwd=2, col = VERMILLION)
legend("bottom",legend = c('Plus','Minus'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n')


plot(conf_beta$minus,ylim=c(15,22),main='Beta experiment',col = BLUE, type = 'b',
     lty = 1, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Beta",
     xlab = "Within phase block of trials")
axis(1, at = 1:4, labels = 1:4)
lines(conf_beta$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2)
error.bar(1:length(conf_beta$minus),conf_beta$minus,conf_beta_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:length(conf_beta$plus),conf_beta$plus,conf_beta_sd$plus,
          lwd=2, col = VERMILLION)
par(mfrow=c(1,1))
dev.off()
Data <- Data_save


