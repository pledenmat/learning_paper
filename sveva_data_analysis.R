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
Data_alpha <- subset(Data,manip=='alpha')
Data_beta <- subset(Data,manip=='beta')

Data$cj <- Data$cj/6
Data$fb <- Data$fb/100
Data$rt2 <- Data$rt + Data$RTconf

setwd(curdir)
# write.csv(Data,file = "alternating_fb_mod.csv",row.names = F)

# test <- read.csv("alternating_fb_mod.csv")

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


# Retrieve fits -----------------------------------------------------------
beta_input <- .1
Nupdate_per_trial <- 1
model <- "allpar"
dt <- .001; sigma <- .1
binning <- F

conditions <- sort(unique(Data$condition))
difflevels <- sort(unique(Data$difflevel)) # Important to sort to match with drift order

ntrial <- Ntrials
par_trace <- data.frame(trial=1:ntrial,sub=rep(subs,each=ntrial*2),
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
  par[par$sub==s,"a0"] <- ldc.results$optim$bestmem[1]
  par[par$sub==s,"b0"] <- ldc.results$optim$bestmem[2]
  par[par$sub==s,"eta"] <- ldc.results$optim$bestmem[4]
  par[par$sub==s,"cost_ldc"] <- ldc.results$optim$bestval

  results <-
    ldc.nn.fit.w(params=c(mean(par[par$sub==s,"a0"]),
                          mean(par[par$sub==s,"b0"]),1,
                          mean(par[par$sub==s,"eta"])),
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


  par[par$sub==s,"manip"] <- unique(temp_dat$manip)
}
summary(par)
par <- par[complete.cases(par$bound),]

# Test DDM fits -----------------------------------------------------------

rm(Simuls_ddm)
nsim <- 1
diffname <- 'difflevel'
difficulty <- difflevels
for (i in 1:Nsub) {
  print(paste('simulating',i,'from',Nsub))
  for (cond in conditions) {
    obs <- subset(Data,sub==subs[i]&condition==cond)
    params <- subset(par,sub==subs[i]&condition==cond)
    # if (!all(complete.cases(params))) {
    #   next
    # }
    drift <- params$drift # Make it more flexible
    params <- as.numeric(as.character(params[1,c("bound","ter","z","vratio")]))
    names(params) <- c('a','ter','z','vratio')
    
    #' First generate DDM predictions
    for (d in 1:length(drift)) {
      if (d == 1) {
        predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(
          v = drift[d], a = params['a'], ter = params['ter'], z = params['z'],
          ntrials = round(nsim*dim(obs)[1]/length(drift)), s = sigma, dt = dt, postdriftmod = params['vratio'], 
          t2distribution = rep(obs[,"RTconf"],times=nsim)
        ))
        names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
        predictions['drift'] <- drift[d]
        predictions[diffname] <- difficulty[d]
      }
      else{
        temp_pred <- data.frame(DDM_with_confidence_slow_fullconfRT(
          v = drift[d], a = params['a'], ter = params['ter'], z = params['z'],
          ntrials = round(nsim*dim(obs)[1]/length(drift)), s = sigma, dt = dt, postdriftmod = params['vratio'], 
          t2distribution = rep(obs[,"RTconf"],times=nsim)
        ))
        names(temp_pred) <- c('rt','resp','cor','evidence2','rt2','cj')
        temp_pred['drift'] <- drift[d]
        temp_pred[diffname] <- difficulty[d]
        predictions <- fastmerge(predictions,temp_pred)
        
      }
    }
    predictions$condition <- cond
    predictions$sub <- subs[i]
    predictions$manip <- unique(obs$manip)
    if (!exists("Simuls_ddm")) {
      Simuls_ddm <- predictions
    } else {
      Simuls_ddm <- rbind(Simuls_ddm,predictions)
    }
  }
}



# Plot trace for the whole experiment -------------------------------------


