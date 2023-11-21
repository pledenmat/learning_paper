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
library(scales)
library(hrbrthemes)
library(ggplot2)
library(geomtextpath)
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

# Extract post-experiment questionnaire
questions <- Data[complete.cases(Data$questionID),
                  c("sub","age","gender","handedness","questionID","questionResp")]

# Remove questionnaire from main data frame
Data <- Data[complete.cases(Data$block),]
rm_col <- c("File","questionID","questionResp","X","ï..sub")
Data <- Data[,-which(names(Data) %in% rm_col)]


subs <- unique(Data$sub); Nsub <- length(subs)
Ntrials <- dim(Data)[1]/Nsub

Data$response <- 0
Data[Data$resp=="['c']",'response'] <- 1

## Diagnostic plot per participant + chance performance testing
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


# Convert RTs to seconds
Data$rt <- Data$rt/1000
Data$RTconf <- Data$RTconf/1000

Data$response[Data$response==0] <- -1

Data$condition <- as.factor(Data$condition)
Data$difflevel <- as.factor(Data$difflevel)
Data$sub <- as.factor(Data$sub)
subs <- unique(Data$sub); Nsub <- length(subs)

Data[Data$manip=='alfa','manip'] <- 'alpha'
names(Data)[1] <- "trial"

# Confidence and feedback on scale from 0 to 1
Data$cj <- Data$cj/6
Data$fb <- Data$fb/100

Data$rt2 <- Data$rt + Data$RTconf # Total RT

# Separate participants according to which condition was first for plotting
Data$group <- "plus_first"
minus_first <- unique(subset(Data,phase==0&condition=='minus')$sub)
Data[Data$sub %in% minus_first,"group"] <- "minus_first"

setwd("..")
Nskip <- 0
Data <- subset(Data,RTconf<5 & rt<5 & rt>.2 & trial>=Nskip)
Ntrials <- Ntrials - Nskip # Because we now skip the first 5 trials
# write.csv(Data,file = "alternating_fb_mod_trim_skip.csv",row.names = F)
write.csv(Data,file = "alternating_fb_mod_trim.csv",row.names = F)
# write.csv(Data,file = "alternating_fb_mod.csv",row.names = F)

Nphase_trial <- length(unique(Data$withinphasetrial))
Nphase_block <- 4
Data$phase_block <-  Data$withinphasetrial %/% (Nphase_trial/Nphase_block)
Data$phase_block <- as.factor(Data$phase_block)
Data$cor <- as.factor(Data$cor)

Data_alpha <- subset(Data,manip=='alpha')
Data_beta <- subset(Data,manip=='beta')

# Demographics
age_alpha <- with(Data_alpha,aggregate(age,by=list(sub),mean))
summary(age_alpha$x)
gender_alpha <- table(Data_alpha$gender)

age_beta <- with(Data_beta,aggregate(age,by=list(sub),mean))
summary(age_beta$x)
gender_beta <- table(Data_beta$gender)

# Behavior analysis -----------------------------------------------------------
if (stat_tests) {
  ## Replication of previous experiments (static effects)
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
  
  ## Evolution of confidence within a switch
  Data_alpha$cor <- as.factor(Data_alpha$cor)
  Data_alpha$phase_block <- as.numeric(Data_alpha$phase_block)
  m.int <- lmer(cj ~ condition*cor*phase_block + (1|sub),
    data=Data_alpha, REML = F, control = lmerControl(optimizer='bobyqa'))
  m.phase <- lmer(cj ~ condition*cor*phase_block + (phase_block|sub),
                data=Data_alpha, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.int,m.phase)
  m.cond <- lmer(cj ~ condition*cor*phase_block + (condition|sub),
                  data=Data_alpha, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.int,m.cond)  
  m.cond.phase <- lmer(cj ~ condition*cor*phase_block + (phase_block+condition|sub),
                 data=Data_alpha, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.cond.phase,m.cond)
  m.condxphase <- lmer(cj ~ condition*cor*phase_block + (phase_block*condition|sub),
                       data=Data_alpha, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.condxphase,m.cond.phase)
  m.condxphase.cor <- lmer(cj ~ condition*cor*phase_block + (phase_block*condition+cor|sub),
                       data=Data_alpha, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.condxphase.cor,m.condxphase)
  m.condxphasexcor <- lmer(cj ~ condition*cor*phase_block + (phase_block*condition+cor+condition:cor|sub),
                           data=Data_alpha, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.condxphasexcor,m.condxphase.cor)
  anova(m.condxphasexcor)
  
  Data_alpha$cor <- as.factor(Data_alpha$cor)
  Data_alpha$phase_block <- as.numeric(Data_alpha$phase_block)
  m.int <- lmer(cj ~ condition*cor*phase_block + (1|sub),
    data=Data_alpha, REML = F, control = lmerControl(optimizer='bobyqa'))
  m.phase <- lmer(cj ~ condition*cor*phase_block + (phase_block|sub),
                data=Data_alpha, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.int,m.phase)
  m.cond <- lmer(cj ~ condition*cor*phase_block + (condition|sub),
                  data=Data_alpha, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.int,m.cond)  
  m.cond.phase <- lmer(cj ~ condition*cor*phase_block + (phase_block+condition|sub),
                 data=Data_alpha, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.cond.phase,m.cond)
  m.condxphase <- lmer(cj ~ condition*cor*phase_block + (phase_block*condition|sub),
                       data=Data_alpha, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.condxphase,m.cond.phase)
  m.condxphase.cor <- lmer(cj ~ condition*cor*phase_block + (phase_block*condition+cor|sub),
                       data=Data_alpha, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.condxphase.cor,m.condxphase)
  m.condxphasexcor <- lmer(cj ~ condition*cor*phase_block + (phase_block*condition+cor+condition:cor|sub),
                           data=Data_alpha, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.condxphasexcor,m.condxphase.cor)
  anova(m.condxphasexcor)
  
  Data_beta$cor <- as.factor(Data_beta$cor)
  Data_beta$phase_block <- as.numeric(Data_beta$phase_block)
  m.beta.int <- lmer(cj ~ condition*cor*phase_block + (1|sub),
                     data=Data_beta, REML = F, control = lmerControl(optimizer='bobyqa'))
  m.beta.phase <- lmer(cj ~ condition*cor*phase_block + (phase_block|sub),
                       data=Data_beta, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.beta.int,m.beta.phase)
  m.beta.cond <- lmer(cj ~ condition*cor*phase_block + (condition|sub),
                      data=Data_beta, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.beta.int,m.beta.cond)  
  m.beta.cond.phase <- lmer(cj ~ condition*cor*phase_block + (phase_block+condition|sub),
                            data=Data_beta, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.beta.cond.phase,m.beta.cond)
  m.beta.condxphase <- lmer(cj ~ condition*cor*phase_block + (phase_block*condition|sub),
                            data=Data_beta, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.beta.condxphase,m.beta.cond.phase)
  m.beta.condxphase.cor <- lmer(cj ~ condition*cor*phase_block + (phase_block*condition+cor|sub),
                                data=Data_beta, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.beta.condxphase.cor,m.beta.condxphase)
  m.beta.condxphasexcor <- lmer(cj ~ condition*cor*phase_block + (phase_block*condition+cor+condition:cor|sub),
                                data=Data_beta, REML = F, control = lmerControl(optimizer='bobyqa'))
  anova(m.beta.condxphasexcor,m.beta.condxphase.cor)
  anova(m.beta.condxphasexcor)
  emm <- emmeans(m.beta.condxphasexcor, ~ phase_block | condition)
  pairs(emm)
  vif(m.beta.condxphasexcor)
  
  

}
test <- Data_alpha[,c("condition","cor","phase_block","cj")]
cor(test)
# Plot Feedback presented ------------------------------------------------
if (plots) {
  go_to("plots")
  go_to("alternating_fb")
  ###
  # Experiment A
  ###
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
  
  jpeg(filename = "feedback_presented.jpg",height = 8,width = 16,units = 'cm', res = 600)
  par(mfrow=c(1,2))
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
  
  ###
  # Experiment B
  ###
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
  dev.off()
  par(mfrow=c(1,1))
}
# Retrieve fits -----------------------------------------------------------
# Some fixed arguments
beta_input <- .1
Nupdate_per_trial <- 1
model <- "allpar"
dt <- .001; sigma <- .1
binning <- F
Nsim <- 20

conditions <- sort(unique(Data$condition))
difflevels <- sort(unique(Data$difflevel)) # Important to sort to match with drift order

totlen <- length(conditions)*length(difflevels)*length(subs)
# DDM parameters + hyperparameters + error

ddm_par <- data.frame(bound = NA, drift = NA, ter = NA,cost_ddm = NA, 
                      sub = rep(subs, each = length(conditions)*length(difflevels)),
                      condition = rep(conditions,each = length(difflevels),
                                      length.out=totlen),
                      difflevel = rep(difflevels, length.out = totlen))

par_trim <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                  sub = rep(subs, each = length(conditions)*length(difflevels)),
                  condition = rep(conditions,each = length(difflevels),
                                  length.out=totlen),
                  difflevel = rep(difflevels, length.out = totlen),
                  manip=NA)
par_trim_skip <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
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
if (!file.exists("anal_sim.Rdata")) {
  anal_sim <- Data[rep(seq_len(nrow(Data)), each=Nsim), c('trial','withinphasetrial','sub','condition','cor')]
  anal_sim$alpha <- NA
  anal_sim$beta <- NA
  anal_sim$sim <- 1:Nsim
}
if (!file.exists("anal_sim_trim.Rdata")) {
  anal_sim_trim <- Data[rep(seq_len(nrow(Data)), each=Nsim), c('trial','withinphasetrial','sub','condition','cor')]
  anal_sim_trim$alpha <- NA
  anal_sim_trim$beta <- NA
  anal_sim_trim$sim <- 1:Nsim
}
if (!file.exists("anal_sim_trim_skip.Rdata")) {
  anal_sim_trim_skip <- Data[rep(seq_len(nrow(Data)), each=Nsim), c('trial','withinphasetrial','sub','condition','cor')]
  anal_sim_trim_skip$alpha <- NA
  anal_sim_trim_skip$beta <- NA
  anal_sim_trim_skip$sim <- 1:Nsim
}
if (!file.exists("anal_sim_beta_learn.Rdata")) {
  anal_sim_beta_learn <- Data[rep(seq_len(nrow(Data)), each=Nsim), c('trial','withinphasetrial','sub','condition','cor')]
  anal_sim_beta_learn$alpha <- NA
  anal_sim_beta_learn$beta <- NA
  anal_sim_beta_learn$cj <- NA
  anal_sim_beta_learn$sim <- 1:Nsim
}
if (!file.exists("anal_sim_alpha_learn.Rdata")) {
  anal_sim_alpha_learn <- Data[rep(seq_len(nrow(Data)), each=Nsim), c('trial','withinphasetrial','sub','condition','cor')]
  anal_sim_alpha_learn$alpha <- NA
  anal_sim_alpha_learn$beta <- NA
  anal_sim_alpha_learn$cj <- NA
  anal_sim_alpha_learn$sim <- 1:Nsim
}
if (!file.exists("anal_sim_both_learn.Rdata")) {
  anal_sim_both_learn <- Data[rep(seq_len(nrow(Data)), each=Nsim), c('trial','withinphasetrial','sub','condition','cor')]
  anal_sim_both_learn$alpha <- NA
  anal_sim_both_learn$beta <- NA
  anal_sim_both_learn$cj <- NA
  anal_sim_both_learn$sim <- 1:Nsim
}
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
  
  ldc_file <- paste0('ldc_nn/trim/batch_',Nupdate_per_trial,'/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
    par_trim[par_trim$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
    par_trim[par_trim$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
    par_trim[par_trim$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
    par_trim[par_trim$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[4]
    par_trim[par_trim$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  }
  
  ldc_file <- paste0('ldc_nn/trim_skip/batch_',Nupdate_per_trial,'/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file)) {
    load(ldc_file)
    par_trim_skip[par_trim_skip$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
    par_trim_skip[par_trim_skip$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
    par_trim_skip[par_trim_skip$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
    par_trim_skip[par_trim_skip$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[4]
    par_trim_skip[par_trim_skip$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  }
  

  ldc_file_beta_learn <- paste0('ldc_nn/trim_skip/beta_learn/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file_beta_learn)) {
    load(ldc_file_beta_learn)
    par_beta_learn[par_beta_learn$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
    par_beta_learn[par_beta_learn$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
    par_beta_learn[par_beta_learn$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
    par_beta_learn[par_beta_learn$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
    par_beta_learn[par_beta_learn$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  }else{
    # # Some fits were not properly saved so retrieved from cluster job output logs
    # ldc_save <- paste0('ldc_nn/trim_skip/beta_learn/slurm-55336494_',s,'.out')
    # if (file.exists(ldc_save)) {
    #   test <- readLines(ldc_save)
    #   check_finish <- "***** summary of DEoptim object ***** "
    #   if (test[length(test)-5] == check_finish) {
    #     fit_par <- unlist(strsplit(test[length(test)-4],split=" "))
    #     a0 <- as.numeric(fit_par[7])
    #     b0 <- as.numeric(fit_par[8])
    #     eta <- as.numeric(fit_par[11])
    #     cost <- as.numeric(unlist(strsplit(test[length(test)-3],split=" "))[8])
    #     ldc.results <- list(optim=list(bestmem=c(a0,b0,1,0,eta),bestval=cost))
    #     save(ldc.results,file=ldc_file_beta_learn)
    #   }
    # }
  }

  ldc_file_alpha_learn <- paste0('ldc_nn/trim_skip/alpha_learn/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file_alpha_learn)) {
    load(ldc_file_alpha_learn)
    par_alpha_learn[par_alpha_learn$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
    par_alpha_learn[par_alpha_learn$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
    par_alpha_learn[par_alpha_learn$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
    par_alpha_learn[par_alpha_learn$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
    par_alpha_learn[par_alpha_learn$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  }
  
  ldc_file_both_learn <- paste0('ldc_nn/trim_skip/both_learn/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file_both_learn)) {
    load(ldc_file_both_learn)
    par_both_learn[par_both_learn$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
    par_both_learn[par_both_learn$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
    par_both_learn[par_both_learn$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
    par_both_learn[par_both_learn$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
    par_both_learn[par_both_learn$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  } else {
    if (s == 23) {
      a0 <- 1.964909
      b0 <- -21.755004
      eta_a <- 0.020717
      eta_b <- 661.525246
      cost <- 0.001528
      ldc.results <- list(optim=list(bestmem=c(a0,b0,1,eta_a,eta_b),bestval=cost))
    }
  }
  
  #' Generate model predictions multiple times to account for the stochastic
  #' nature of estimating single trial accumulated evidence 
  if (file.exists("anal_sim_trim.Rdata")) {
    if (!exists("anal_sim_trim")) {
      load("anal_sim_trim.Rdata")
    }
  } else {
    for (i in 1:Nsim) {
      results <-
        ldc.nn.fit.w(params=c(mean(par_trim[par_trim$sub==subs[s],"a0"]),
                                   mean(par_trim[par_trim$sub==subs[s],"b0"]),1,
                                   mean(par_trim[par_trim$sub==subs[s],"eta_a"])),
                     ddm_params = ddm_params,
                     obs=temp_dat,returnFit = F,
                     Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                     dt = dt, sigma = sigma)
      
      anal_sim_trim[anal_sim_trim$sim==i&anal_sim_trim$sub==subs[s] ,'cj'] <- results$pred
      anal_sim_trim[anal_sim_trim$sim==i&anal_sim_trim$sub==subs[s] ,'alpha'] <- results$trace[,1]
      anal_sim_trim[anal_sim_trim$sim==i&anal_sim_trim$sub==subs[s] ,'beta'] <- results$trace[,2]  
    }
  }
  
  if (file.exists("anal_sim_trim_skip.Rdata")) {
    if (!exists("anal_sim_trim_skip")) {
      load("anal_sim_trim_skip.Rdata")
    }
  } else {
    for (i in 1:Nsim) {
      results <-
        ldc.nn.fit.w(params=c(mean(par_trim_skip[par_trim_skip$sub==subs[s],"a0"]),
                                        mean(par_trim_skip[par_trim_skip$sub==subs[s],"b0"]),1,
                                        mean(par_trim_skip[par_trim_skip$sub==subs[s],"eta_a"])),
                     ddm_params = ddm_params,
                     obs=temp_dat,returnFit = F,
                     Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                     dt = dt, sigma = sigma)
      
      anal_sim_trim_skip[anal_sim_trim_skip$sim==i&anal_sim_trim_skip$sub==subs[s] ,'cj'] <- results$pred
      anal_sim_trim_skip[anal_sim_trim_skip$sim==i&anal_sim_trim_skip$sub==subs[s] ,'alpha'] <- results$trace[,1]
      anal_sim_trim_skip[anal_sim_trim_skip$sim==i&anal_sim_trim_skip$sub==subs[s] ,'beta'] <- results$trace[,2]  
    }
  }
  
  #' Do the same with the alternative models
  if (file.exists("anal_sim_beta_learn.Rdata")) {
    if (!exists("anal_sim_beta_learn")) {
      load("anal_sim_beta_learn.Rdata")
    }
  } else {
    for (i in 1:Nsim) {
      results <-
        ldc.nn.fit.w(params=c(mean(par_beta_learn[par_beta_learn$sub==subs[s],"a0"]),
                              mean(par_beta_learn[par_beta_learn$sub==subs[s],"b0"]),1,0,
                              mean(par_beta_learn[par_beta_learn$sub==subs[s],"eta_b"])),
                     ddm_params = ddm_params,
                     obs=temp_dat,returnFit = F,eta_sep=T,
                     Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                     dt = dt, sigma = sigma)
      anal_sim_beta_learn[anal_sim_beta_learn$sim==i&anal_sim_beta_learn$sub==subs[s] ,'cj'] <- results$pred
      anal_sim_beta_learn[anal_sim_beta_learn$sim==i&anal_sim_beta_learn$sub==subs[s] ,'alpha'] <- results$trace[,1]
      anal_sim_beta_learn[anal_sim_beta_learn$sim==i&anal_sim_beta_learn$sub==subs[s] ,'beta'] <- results$trace[,2]  
    }
  }
  
  if (file.exists("anal_sim_alpha_learn.Rdata")) {
    load("anal_sim_alpha_learn.Rdata")
  } else {
    for (i in 1:Nsim) {
      results <-
        ldc.nn.fit.w(params=c(mean(par_alpha_learn[par_alpha_learn$sub==subs[s],"a0"]),
                              mean(par_alpha_learn[par_alpha_learn$sub==subs[s],"b0"]),1,
                              mean(par_alpha_learn[par_alpha_learn$sub==subs[s],"eta_a"]),0),
                     ddm_params = ddm_params,
                     obs=temp_dat,returnFit = F,eta_sep=T,
                     Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                     dt = dt, sigma = sigma)
      
      anal_sim_alpha_learn[anal_sim_alpha_learn$sim==i&anal_sim_alpha_learn$sub==subs[s] ,'cj'] <- results$pred
      anal_sim_alpha_learn[anal_sim_alpha_learn$sim==i&anal_sim_alpha_learn$sub==subs[s] ,'alpha'] <- results$trace[,1]
      anal_sim_alpha_learn[anal_sim_alpha_learn$sim==i&anal_sim_alpha_learn$sub==subs[s] ,'beta'] <- results$trace[,2]  
    }
  }

  if (file.exists("anal_sim_both_learn.Rdata")) {
    if (!exists("anal_sim_both_learn")) {
      load("anal_sim_both_learn.Rdata")
    }
  } else {
    for (i in 1:Nsim) {
      results <-
        ldc.nn.fit.w(params=c(mean(par_both_learn[par_both_learn$sub==subs[s],"a0"]),
                              mean(par_both_learn[par_both_learn$sub==subs[s],"b0"]),1,
                              mean(par_both_learn[par_both_learn$sub==subs[s],"eta_a"]),
                                   mean(par_both_learn[par_both_learn$sub==subs[s],"eta_b"])),
                              ddm_params = ddm_params,
                              obs=temp_dat,returnFit = F,eta_sep=T,
                              Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                              dt = dt, sigma = sigma)
                     anal_sim_both_learn[anal_sim_both_learn$sim==i&anal_sim_both_learn$sub==subs[s] ,'cj'] <- results$pred
                     anal_sim_both_learn[anal_sim_both_learn$sim==i&anal_sim_both_learn$sub==subs[s] ,'alpha'] <- results$trace[,1]
                     anal_sim_both_learn[anal_sim_both_learn$sim==i&anal_sim_both_learn$sub==subs[s] ,'beta'] <- results$trace[,2]  
    }
  }

  par_trim[par_trim$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  par_trim_skip[par_trim_skip$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  par_both_learn[par_both_learn$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  par_beta_learn[par_beta_learn$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  par_alpha_learn[par_alpha_learn$sub==subs[s],"manip"] <- unique(temp_dat$manip)
}
# Save model simulations
if (!file.exists("anal_sim.Rdata")) {
  save(anal_sim,file="anal_sim.Rdata")
}
if (!file.exists("anal_sim_trim.Rdata")) {
  save(anal_sim_trim,file="anal_sim_trim.Rdata")
}
if (!file.exists("anal_sim_trim_skip.Rdata")) {
  save(anal_sim_trim_skip,file="anal_sim_trim_skip.Rdata")
}
if (!file.exists("anal_sim_beta_learn.Rdata")) {
  save(anal_sim_beta_learn,file="anal_sim_beta_learn.Rdata")
}
if (!file.exists("anal_sim_alpha_learn.Rdata")) {
  save(anal_sim_alpha_learn,file="anal_sim_alpha_learn.Rdata")
}
if (!file.exists("anal_sim_both_learn.Rdata")) {
  save(anal_sim_both_learn,file="anal_sim_both_learn.Rdata")
}

# Merge parameters from every model
par_trim$model <- "trim"
par_trim_skip$model <- "trim_skip"
par_alpha_learn$model <- "alpha"
par_beta_learn$model <- "beta"
par_both_learn$model <- "both"

par <- rbind(par_alpha_learn,par_beta_learn,par_both_learn,par_trim,par_trim_skip)
par <- merge(par,ddm_par,by = c("sub","condition","difflevel"))

# Add confidence prediction and parameter trace from each model to the empirical data frame
cj_pred <- with(anal_sim,aggregate(cj,by=list(trial,sub),mean))
alpha <- with(anal_sim,aggregate(alpha,by=list(trial,sub),mean))
beta <- with(anal_sim,aggregate(beta,by=list(trial,sub),mean))
names(cj_pred) <- c("trial","sub","cj_pred")
names(alpha) <- c("trial","sub","alpha")
names(beta) <- c("trial","sub","beta")
anal_sim_mean <- merge(merge(cj_pred,alpha),beta)
Data <- Data[, names(Data) != "cj_pred"]
Data <- merge(Data,anal_sim_mean)

cj_pred <- with(anal_sim_trim,aggregate(cj,by=list(trial,sub),mean))
alpha <- with(anal_sim_trim,aggregate(alpha,by=list(trial,sub),mean))
beta <- with(anal_sim_trim,aggregate(beta,by=list(trial,sub),mean))
names(cj_pred) <- c("trial","sub","cj_pred_trim")
names(alpha) <- c("trial","sub","alpha_trim")
names(beta) <- c("trial","sub","beta_trim")
anal_sim_trim_mean <- merge(merge(cj_pred,alpha),beta)
Data <- merge(Data,anal_sim_trim_mean)

cj_pred <- with(anal_sim_trim_skip,aggregate(cj,by=list(trial,sub),mean))
alpha <- with(anal_sim_trim_skip,aggregate(alpha,by=list(trial,sub),mean))
beta <- with(anal_sim_trim_skip,aggregate(beta,by=list(trial,sub),mean))
names(cj_pred) <- c("trial","sub","cj_pred_trim_skip")
names(alpha) <- c("trial","sub","alpha_trim_skip")
names(beta) <- c("trial","sub","beta_trim_skip")
anal_sim_trim_skip_mean <- merge(merge(cj_pred,alpha),beta)
Data <- merge(Data,anal_sim_trim_skip_mean)

cj_pred <- with(anal_sim_beta_learn,aggregate(cj,by=list(trial,sub),mean))
alpha <- with(anal_sim_beta_learn,aggregate(alpha,by=list(trial,sub),mean))
beta <- with(anal_sim_beta_learn,aggregate(beta,by=list(trial,sub),mean))
names(cj_pred) <- c("trial","sub","cj_pred_beta_learn")
names(alpha) <- c("trial","sub","alpha_beta_learn")
names(beta) <- c("trial","sub","beta_beta_learn")
anal_sim_beta_learn_mean <- merge(merge(cj_pred,alpha),beta)
Data <- merge(Data,anal_sim_beta_learn_mean)


cj_pred <- with(anal_sim_alpha_learn,aggregate(cj,by=list(trial,sub),mean))
alpha <- with(anal_sim_alpha_learn,aggregate(alpha,by=list(trial,sub),mean))
beta <- with(anal_sim_alpha_learn,aggregate(beta,by=list(trial,sub),mean))
names(cj_pred) <- c("trial","sub","cj_pred_alpha_learn")
names(alpha) <- c("trial","sub","alpha_alpha_learn")
names(beta) <- c("trial","sub","beta_alpha_learn")
anal_sim_alpha_learn_mean <- merge(merge(cj_pred,alpha),beta)
Data <- merge(Data,anal_sim_alpha_learn_mean)

cj_pred <- with(anal_sim_both_learn,aggregate(cj,by=list(trial,sub),mean))
alpha <- with(anal_sim_both_learn,aggregate(alpha,by=list(trial,sub),mean))
beta <- with(anal_sim_both_learn,aggregate(beta,by=list(trial,sub),mean))
names(cj_pred) <- c("trial","sub","cj_pred_both_learn")
names(alpha) <- c("trial","sub","alpha_both_learn")
names(beta) <- c("trial","sub","beta_both_learn")
anal_sim_both_learn_mean <- merge(merge(cj_pred,alpha),beta)
Data <- merge(Data,anal_sim_both_learn_mean)

Data_alpha <- subset(Data,manip=='alpha')
Data_beta <- subset(Data,manip=='beta')

# Compute rolling mean per subject ----------------------------------------
n <- 25 # Rolling mean window size
n_err <- 25

trial_conf_sub <- with(Data,aggregate(cj,by=list(trial,cor,sub),mean))
names(trial_conf_sub) <- c("trial","cor","sub","cj")

pred_conf_sub_beta_learn <- with(Data,aggregate(cj_pred_beta_learn,by=list(trial,cor,sub),mean))
names(pred_conf_sub_beta_learn) <- c("trial","cor","sub","cj")

pred_conf_sub_alpha_learn <- with(Data,aggregate(cj_pred_alpha_learn,by=list(trial,cor,sub),mean))
names(pred_conf_sub_alpha_learn) <- c("trial","cor","sub","cj")

pred_conf_sub_both_learn <- with(Data,aggregate(cj_pred_both_learn,by=list(trial,cor,sub),mean))
names(pred_conf_sub_both_learn) <- c("trial","cor","sub","cj")

pred_conf_sub_trim_skip <- with(Data,aggregate(cj_pred_trim_skip,by=list(trial,cor,sub),mean))
names(pred_conf_sub_trim_skip) <- c("trial","cor","sub","cj")

pred_conf_sub_trim <- with(Data,aggregate(cj_pred_trim,by=list(trial,cor,sub),mean))
names(pred_conf_sub_trim) <- c("trial","cor","sub","cj")

trials <- data.frame(trial=rep((0:(Ntrials-1))+Nskip,each=2),
                     cor=c(0,1),sub=rep(subs,each=Ntrials*2))

cj_ma <- merge(trial_conf_sub,trials,all=T)
cj_pred_ma_beta_learn <- merge(pred_conf_sub_beta_learn,trials,all=T)
cj_pred_ma_alpha_learn <- merge(pred_conf_sub_alpha_learn,trials,all=T)
cj_pred_ma_both_learn <- merge(pred_conf_sub_both_learn,trials,all=T)
cj_pred_ma_trim <- merge(pred_conf_sub_trim,trials,all=T)
cj_pred_ma_trim_skip <- merge(pred_conf_sub_trim_skip,trials,all=T)

ma <- function(x,n,names){
  return(rollapply(x[,names], width=n, FUN=function(x) mean(x, na.rm=TRUE),partial=TRUE, align="center"))
}
for (s in subs) {
  cj_ma[cj_ma$sub==s&cj_ma$cor==0,"cj"] <- ma(subset(cj_ma,sub==s&cor==0),n_err,"cj")
  cj_ma[cj_ma$sub==s&cj_ma$cor==1,"cj"] <- ma(subset(cj_ma,sub==s&cor==1),n,"cj")
  cj_pred_ma_beta_learn[cj_pred_ma_beta_learn$sub==s&cj_pred_ma_beta_learn$cor==0,"cj"] <- ma(subset(cj_pred_ma_beta_learn,sub==s&cor==0),n_err,"cj")
  cj_pred_ma_beta_learn[cj_pred_ma_beta_learn$sub==s&cj_pred_ma_beta_learn$cor==1,"cj"] <- ma(subset(cj_pred_ma_beta_learn,sub==s&cor==1),n,"cj")
  cj_pred_ma_alpha_learn[cj_pred_ma_alpha_learn$sub==s&cj_pred_ma_alpha_learn$cor==0,"cj"] <- ma(subset(cj_pred_ma_alpha_learn,sub==s&cor==0),n_err,"cj")
  cj_pred_ma_alpha_learn[cj_pred_ma_alpha_learn$sub==s&cj_pred_ma_alpha_learn$cor==1,"cj"] <- ma(subset(cj_pred_ma_alpha_learn,sub==s&cor==1),n,"cj")
  cj_pred_ma_both_learn[cj_pred_ma_both_learn$sub==s&cj_pred_ma_both_learn$cor==0,"cj"] <- ma(subset(cj_pred_ma_both_learn,sub==s&cor==0),n_err,"cj")
  cj_pred_ma_both_learn[cj_pred_ma_both_learn$sub==s&cj_pred_ma_both_learn$cor==1,"cj"] <- ma(subset(cj_pred_ma_both_learn,sub==s&cor==1),n,"cj")
  cj_pred_ma_trim_skip[cj_pred_ma_trim_skip$sub==s&cj_pred_ma_trim_skip$cor==0,"cj"] <- ma(subset(cj_pred_ma_trim_skip,sub==s&cor==0),n_err,"cj")
  cj_pred_ma_trim_skip[cj_pred_ma_trim_skip$sub==s&cj_pred_ma_trim_skip$cor==1,"cj"] <- ma(subset(cj_pred_ma_trim_skip,sub==s&cor==1),n,"cj")
  cj_pred_ma_trim[cj_pred_ma_trim$sub==s&cj_pred_ma_trim$cor==0,"cj"] <- ma(subset(cj_pred_ma_trim,sub==s&cor==0),n_err,"cj")
  cj_pred_ma_trim[cj_pred_ma_trim$sub==s&cj_pred_ma_trim$cor==1,"cj"] <- ma(subset(cj_pred_ma_trim,sub==s&cor==1),n,"cj")
}
# Model comparison --------------------------------------------------------
### Residuals
model_residuals <- data.frame(sub=par$sub,manip=par$manip,
                              cost_baseline=par$cost_ldc,
                              cost_alpha=par_alpha_learn$cost_ldc,
                              cost_beta=par_beta_learn$cost_ldc)
table(apply(model_residuals[,3:5],1,which.min))
jpeg(filename="comparison_residuals.jpg",width=10,height=8,units='cm',res=300)
plot(density(model_residuals$cost_beta,from = 0,to = .02), bty = 'n', main = "", xlab = "Residuals",
     col = "darkred",lwd=2)
lines(density(model_residuals$cost_baseline,from = 0,to = .02),col="orange",lwd=2)
lines(density(model_residuals$cost_alpha,from = 0,to = .02),col="cyan",lwd=2)
legend("topright",c("Baseline","Alpha learn","Beta learn"),col=c("orange","cyan","darkred"),
       lty = 1,bty='n',lwd=2)
dev.off()

### Correlations with empirical
cor_baseline <- with(Data,aggregate(cj,by=list(sub,manip),mean))
names(cor_baseline) <- c('sub','manip','cor')
cor_baseline$cor <- NA
cor_alpha_learn <- cor_baseline
names(cor_alpha_learn) <- c('sub','manip','cor_alpha')
cor_beta_learn <- cor_baseline
names(cor_beta_learn) <- c('sub','manip','cor_beta')
for (s in 1:Nsub) {
  tempdat <- subset(Data,sub==subs[s])
  cor_baseline[cor_baseline$sub==subs[s],'cor'] <- cor(tempdat$cj,tempdat$cj_pred)
  cor_alpha_learn[cor_alpha_learn$sub==subs[s],'cor_alpha'] <- cor(tempdat$cj,tempdat$cj_pred_alpha_learn)
  cor_beta_learn[cor_beta_learn$sub==subs[s],'cor_beta'] <- cor(tempdat$cj,tempdat$cj_pred_beta_learn)
}
cors <- merge(merge(cor_baseline,cor_alpha_learn),cor_beta_learn)
table(apply(cors[,3:5],1,which.min))

jpeg(filename="comparison_correlation.jpg",width=10,height=8,units='cm',res=300)
plot(density(cors$cor_beta,from = 0,to = 1), bty = 'n', main = "", xlab = "Empirical/Prediction correlation",
     col = "darkred",lwd=2)
lines(density(cors$cor,from = 0,to = 1),col="orange",lwd=2)
lines(density(cors$cor_alpha,from = 0,to = 1),col="cyan",lwd=2)
legend("topleft",c("Baseline","Alpha learn","Beta learn"),col=c("orange","cyan","darkred"),
       lty = 1,bty='n',lwd=2)
dev.off()
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
if (stat_tests) {
  
  m.int <- lmer(data = Data_alpha, cj ~ condition*phase_block*cor + (1|sub),REML = F)
  m.cond <- lmer(data = Data_alpha, cj ~ condition*phase_block*cor + (condition|sub),REML = F)
  anova(m.int,m.cond)
  m.cond.cor <- lmer(data = Data_alpha, cj ~ condition*phase_block*cor + (condition+cor|sub),REML = F)
  anova(m.cond,m.cond.cor)
  anova(m.cond.cor)
  emm <- emmeans(m.cond.cor, ~ condition | phase_block)
  pairs(emm)
  emm <- emmeans(m.cond.cor, ~ phase_block | condition)
  pairs(emm)
  Data_alpha$phase_block <- as.numeric(Data_alpha$phase_block)
  m.cond.cor.cont <- lmer(data = Data_alpha, cj ~ condition*phase_block*cor + (condition+cor|sub),REML = F)
  anova(m.cond.cor.cont)
  emm <- emmeans(m.cond.cor.cont, ~ condition | phase_block)
  pairs(emm)
  emm <- emmeans(m.cond.cor.cont, ~ condition)
  pairs(emm)
  emm <- emmeans(m.cond.cor.cont, ~ phase_block)
  pairs(emm)
  m.cond.cor.cont2 <- lmer(data = subset(Data_alpha,phase>0), cj ~ condition*phase_block*cor + (condition+cor|sub),REML = F)
  anova(m.cond.cor.cont2)
  
  m.int <- lmer(data = Data_beta, cj ~ condition*phase_block*cor + (1|sub),REML = F)
  m.cond <- lmer(data = Data_beta, cj ~ condition*phase_block*cor + (condition|sub),REML = F)
  anova(m.int,m.cond)
  m.cond.cor <- lmer(data = Data_beta, cj ~ condition*phase_block*cor + (condition+cor|sub),REML = F)
  anova(m.cond,m.cond.cor)
  anova(m.cond.cor)
  emm <- emmeans(m.cond.cor, ~ condition | cor | phase_block)
  pairs(emm)
}

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
tiff('trace_aggreg_cor.tiff', width = 36, height = 9, units = 'cm', res = 300)
layout(matrix(c(1,3,2,4),ncol=2),heights=c(1,7))
par(mar=c(0,0,0,0))
plot.new()
title(cex.main=cex.title*.66/.83,line=title_line+1,main = expression(paste(alpha,"-Manipulated Feedback")))
plot.new()
title(cex.main=cex.title*.66/.83,line=title_line+1,main = expression(paste(beta,"-Manipulated Feedback")))
par(mar=c(4,5,0,1.1)+.1)

plot(conf_alpha$minus_0,ylim=c(.7,.9),main='',col = BLUE, type = 'b',
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",
     xlab = "Trials",cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:4, labels = c("[1:32]","[33:63]","[64:94]","[95:126]"),cex.axis=cex.axis*.66/.83)
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
legend("bottom",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)


plot(conf_beta$minus_0,ylim=c(.7,.9),main='',col = BLUE, type = 'b',
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",
     xlab = "Trials",cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:4, labels = c("[1:32]","[33:63]","[64:94]","[95:126]"),cex.axis=cex.axis*.66/.83)
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
if (stat_tests) {
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
  
  m.alpha.b <- lmer(data = Data_alpha, beta_beta_learn ~ phase_block*condition + (condition + phase_block|sub), REML = F,
                    control = lmerControl(optimizer = "bobyqa"))
  anova(m.alpha.b)
  emm <- emmeans(m.alpha.b, ~ condition)
  pairs(emm)
  
  # Do the same with beta experiment
  m.beta.a <- lmer(data = Data_beta, alpha ~ phase_block*condition + (phase_block+condition|sub), REML = F,
                   control = lmerControl(optimizer = "bobyqa"))
  anova(m.beta.a)
  emm <- emmeans(m.beta.a, ~ condition | phase_block)
  pairs(emm)
  emm <- emmeans(m.beta.a, ~ phase_block | condition)
  pairs(emm)
  
  m.beta.b <- lmer(data = Data_beta, beta ~ phase_block*condition + (condition + phase_block|sub), REML = F,
                   control = lmerControl(optimizer = "bobyqa"))
  anova(m.beta.b)
  emm <- emmeans(m.beta.b, ~ condition | phase_block)
  pairs(emm)
  emm <- emmeans(m.beta.b, ~ phase_block | condition)
  pairs(emm)
  
  m.beta.b <- lmer(data = Data_beta, beta_beta_learn ~ phase_block*condition + (condition + phase_block|sub), REML = F,
                   control = lmerControl(optimizer = "bobyqa"))
  anova(m.beta.b)
  emm <- emmeans(m.beta.b, ~ condition | phase_block)
  pairs(emm)
  emm <- emmeans(m.beta.b, ~ phase_block | condition)
  pairs(emm)
  
  Data_alpha$phase_block <- as.numeric(Data_alpha$phase_block)
  Data_beta$phase_block <- as.numeric(Data_beta$phase_block)
  
  m.alpha.b.cont <- lmer(data = Data_alpha, beta_beta_learn ~ phase_block*condition + (condition + phase_block|sub), REML = F,
                    control = lmerControl(optimizer = "bobyqa"))
  anova(m.alpha.b.cont)
  anova(m.alpha.b)
  m.beta.b.cont <- lmer(data = Data_beta, beta_beta_learn ~ phase_block*condition + (condition + phase_block|sub), REML = F,
                         control = lmerControl(optimizer = "bobyqa"))
  anova(m.beta.b.cont)
  emm <- emmeans(m.beta.b.cont, ~ phase_block * condition)
  pairs(emm)
  contrasts <- contrast(emm, interaction = "pairwise")
  summary(contrasts)
  summary(m.beta.b)

}

# Plot all this
tiff('param_trace_aggreg.tiff', width = 36, height = 14, units = 'cm', res = 300)
layout(matrix(c(1,3,5,2,4,6),ncol=2),heights = c(1,6.5,6.5))

par(mar=c(0,0,0,0))
plot.new()
title(cex.main=cex.title,line=title_line,main = expression(paste(alpha,"-Manipulated Feedback")))
plot.new()
title(cex.main=cex.title,line=title_line,main = expression(paste(beta,"-Manipulated Feedback")))
par(mar=c(4,5,0,1.1)+.1)

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

plot(conf_alpha$minus,ylim=c(5,13),main=NULL,col = BLUE, type = 'b',
     lty = 1, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Alpha",
     xlab = "Trials",cex.lab=cex.lab,cex.axis=cex.axis)
axis(1, at = 1:4, labels = c("[1:32]","[33:63]","[64:94]","[95:126]"),cex.axis=cex.axis)
lines(conf_alpha$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2)
error.bar(1:length(conf_alpha$minus),conf_alpha$minus,conf_alpha_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:length(conf_alpha$plus),conf_alpha$plus,conf_alpha_sd$plus,
          lwd=2, col = VERMILLION)
legend("bottom",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T,cex=cex.legend, bty = 'n')


plot(conf_beta$minus,ylim=c(5,13),main=NULL,col = BLUE, type = 'b',
     lty = 1, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Alpha",
     xlab = "Trials",cex.lab=cex.lab,cex.axis=cex.axis)
axis(1, at = 1:4, labels = c("[1:32]","[33:63]","[64:94]","[95:126]"),cex.axis=cex.axis)
lines(conf_beta$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2)
error.bar(1:length(conf_beta$minus),conf_beta$minus,conf_beta_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:length(conf_beta$plus),conf_beta$plus,conf_beta_sd$plus,
          lwd=2, col = VERMILLION)

Nalpha <- length(unique(Data_alpha$sub))
Nbeta <- length(unique(Data_beta$sub))
conf_group <- with(Data,aggregate(beta,by=list(phase_block,manip,condition),mean))
conf_group_sd <- with(Data,aggregate(beta,by=list(phase_block,manip,condition),sd))
names(conf_group) <- c('phase_block','manip','condition','cj')
names(conf_group_sd) <- c('phase_block','manip','condition','cj')
conf_group$cj <- conf_group$cj/10
conf_group_sd$cj <- conf_group_sd$cj/10

conf_alpha <- cast(subset(conf_group,manip=='alpha'),phase_block~condition)
conf_beta <- cast(subset(conf_group,manip=='beta'),phase_block~condition)
conf_alpha_sd <- subset(conf_group_sd,manip=='alpha')
conf_alpha_sd$cj <- conf_alpha_sd$cj/sqrt(Nalpha)
conf_alpha_sd <- cast(conf_alpha_sd, phase_block ~ condition)
conf_beta_sd <- subset(conf_group_sd,manip=='beta')
conf_beta_sd$cj <- conf_beta_sd$cj/sqrt(Nbeta)
conf_beta_sd <- cast(conf_beta_sd, phase_block ~ condition)

plot(conf_alpha$minus,ylim=c(1.5,2.2),main=NULL,col = BLUE, type = 'b',
     lty = 1, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Beta",
     xlab = "Trials",cex.axis=cex.axis,cex.lab=cex.lab)
axis(1, at = 1:4, labels = c("[1:32]","[33:63]","[64:94]","[95:126]"),cex.axis=cex.axis)
lines(conf_alpha$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2)
error.bar(1:length(conf_alpha$minus),conf_alpha$minus,conf_alpha_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:length(conf_alpha$plus),conf_alpha$plus,conf_alpha_sd$plus,
          lwd=2, col = VERMILLION)
legend("bottom",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T,cex=cex.legend, bty = 'n')


plot(conf_beta$minus,ylim=c(1.5,2.2),main=NULL,col = BLUE, type = 'b',
     lty = 1, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Beta",
     xlab = "Trials",cex.axis=cex.axis,cex.lab=cex.lab)
axis(1, at = 1:4, labels = c("[1:32]","[33:63]","[64:94]","[95:126]"),cex.axis=cex.axis)
lines(conf_beta$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2)
error.bar(1:length(conf_beta$minus),conf_beta$minus,conf_beta_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:length(conf_beta$plus),conf_beta$plus,conf_beta_sd$plus,
          lwd=2, col = VERMILLION)
par(mfrow=c(1,1))
dev.off()
# Not in Manuscript but cool ----------------------------------------------
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



# Plot Accuracy ~ Confidence ----------------------------------------------

cor_cj <- with(Data,aggregate(cor,by=list(cj,sub),mean))
names(cor_cj) <- c('cj','sub','cor')
cor_cj$cj <- cor_cj$cj*6
count_data <- table(cor_cj$cj)
r <- cor(cor_cj$cj,cor_cj$cor)
cor_cj <- as.matrix(cast(cor_cj, sub~cj))
cexlab <- 2.25
cexax <- 1.5
cexmain <- 2.25
cexpt <- 1.5
ntick <- 2
tiff("accuracy_confidence.tiff",width = 10.5, height = 10.5,units='cm',res=300)
par(mar=c(4,4,4,1)+.1)
plot(colMeans(cor_cj,na.rm=T),ylim=c(0,1),xlab='',col='white',
     ylab='',bty='n', main = paste0('r = ',round(r,3)),
     cex.axis=cexax,cex.main=cexmain,yaxt='n')
mtext(side=1,line=2.5,text='Confidence',cex=cexlab)
mtext(side=2,line=2.5,text='Accuracy',cex=cexlab)
axis(2,at=0:ntick/ntick,labels = 0:ntick/ntick,cex.axis = cexax)
for (i in 1:6) {
  points(jitter(rep(i,nrow(cor_cj)),amount=.1),cor_cj[,i],pch=16,col=rgb(.5,.5,.5,.2),cex=cexpt)
}
lines(colMeans(cor_cj,na.rm=T),type='b',pch=16,lwd=3,cex=cexpt)
error.bar(1:6,colMeans(cor_cj,na.rm=T),colSds(cor_cj,na.rm=T)/sqrt(count_data),lwd=3)
dev.off()
par(mar=c(5,4,4,2)+.1)


# Analyze variability in parameter trace btw simulations ------------------
go_to("param_trace_per_sub")
col_indiv <- rgb(.7,.7,.7,.2)
for (s in subs) {
  print(paste("Plotting parameter trace of sub",s))
  temp_alpha <- as.matrix(cast(subset(anal_sim,sub==s),trial~sim, value = 'alpha'))
  temp_beta <- as.matrix(cast(subset(anal_sim,sub==s),trial~sim, value = 'beta'))
  temp_cj <- as.matrix(cast(subset(anal_sim,sub==s),trial~sim, value = 'cj'))
  jpeg(filename=paste0("param_trace_",s,".jpg"),width=45,height=15,units='cm',res=300)
  par(mfrow=c(1,3))
  plot(rowMeans(temp_alpha,na.rm=T),col='white',xlab='Trial',ylab='Alpha',
       bty='n',ylim=range(temp_alpha,na.rm=T))
  for (i in 1:Nsim) {
    lines(temp_alpha[,i], col = col_indiv)
  }
  lines(rowMeans(temp_alpha,na.rm=T),type='l',lwd=1)
  plot(rowMeans(temp_beta,na.rm=T),col='white',xlab='Trial',ylab='beta',
       bty='n',ylim=range(temp_beta,na.rm=T))
  for (i in 1:Nsim) {
    lines(temp_beta[,i], col = col_indiv)
  }
  lines(rowMeans(temp_beta,na.rm=T),type='l',lwd=1)
  plot(rowMeans(temp_cj,na.rm=T),col='white',xlab='Trial',ylab='Confidence',
       bty='n',ylim=range(temp_cj,na.rm=T))
  for (i in 1:Nsim) {
    lines(temp_cj[,i], col = col_indiv)
  }
  lines(rowMeans(temp_cj,na.rm=T),type='l',lwd=1)
  
  dev.off()
}
setwd("..")
par(mfrow=c(1,1))



# Trash -------------------------------------------------------------------
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


# Plot traces trim skip ---------------------------------------------------

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

jpeg(filename = "traces_trim_skip.jpg",units = 'cm',width = 42,height = 30,res=300)
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

# Plot ExpA trace
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


# Plot ExpA predictions
plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)

conf_min <- with(subset(cj_pred_ma_trim_skip,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_trim_skip,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_trim_skip,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_trim_skip,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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

# Plot ExpB predictions 
plus_first <- unique(subset(Data_beta,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_beta,phase==0&condition=='minus')$sub)

conf_min <- with(subset(cj_pred_ma_trim_skip,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_trim_skip,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_trim_skip,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_trim_skip,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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

# Plot parameter traces 
par_higheta <- subset(par,eta>10000)

plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)

# Plot trace alpha experiment
alpha_trace <- with(Data_alpha,aggregate(alpha_trim_skip,by=list(sub=sub,trial=trial,condition=condition),mean))
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
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_minus,na.rm=T) + 
                                   colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(alpha_trace_minus,na.rm=T) - 
                                                                                          colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_plus,na.rm=T) + 
                                   colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(alpha_trace_plus,na.rm=T) - 
                                                                                        colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))


beta_trace <- with(Data_alpha,aggregate(beta_trim_skip,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(0,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(beta_trace_plus,na.rm = T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_minus,na.rm=T) + 
                                   colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(beta_trace_minus,na.rm=T) - 
                                                                                         colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_plus,na.rm=T) + 
                                   colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(beta_trace_plus,na.rm=T) - 
                                                                                       colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))

# Plot trace Beta experiment

plus_first <- unique(subset(Data_beta,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_beta,phase==0&condition=='minus')$sub)

alpha_trace <- with(Data_beta,aggregate(alpha_trim_skip,by=list(sub=sub,trial=trial,condition=condition),mean))
alpha_trace_minus <- cast(subset(alpha_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
count_plus <- count_plus[2:length(count_plus)]
count_minus <- count_minus[2:length(count_minus)]

plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
     ylim=c(-10,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(alpha_trace_plus,na.rm=T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_minus,na.rm=T) + 
                                   colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(alpha_trace_minus,na.rm=T) - 
                                                                                          colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_plus,na.rm=T) + 
                                   colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(alpha_trace_plus,na.rm=T) - 
                                                                                        colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))


beta_trace <- with(Data_beta,aggregate(beta_trim_skip,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(0,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(beta_trace_plus,na.rm = T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_minus,na.rm=T) + 
                                   colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(beta_trace_minus,na.rm=T) - 
                                                                                         colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_plus,na.rm=T) + 
                                   colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(beta_trace_plus,na.rm=T) - 
                                                                                       colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))

dev.off()


# Plot traces trim model --------------------------------------------------
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

jpeg(filename = "traces_trim.jpg",units = 'cm',width = 42,height = 30,res=300)
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

## Plot ExpA trace

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


# Plot ExpA predictions 
plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)

conf_min <- with(subset(cj_pred_ma_trim,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_trim,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_trim,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_trim,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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

# Plot ExpB predictions 
plus_first <- unique(subset(Data_beta,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_beta,phase==0&condition=='minus')$sub)

conf_min <- with(subset(cj_pred_ma_trim,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_trim,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_trim,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_trim,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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

# Plot parameter traces 
par_higheta <- subset(par,eta>10000)

plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)

# Plot trace alpha experiment
alpha_trace <- with(Data_alpha,aggregate(alpha_trim,by=list(sub=sub,trial=trial,condition=condition),mean))
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
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_minus,na.rm=T) + 
                                   colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(alpha_trace_minus,na.rm=T) - 
                                                                                          colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_plus,na.rm=T) + 
                                   colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(alpha_trace_plus,na.rm=T) - 
                                                                                        colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))


beta_trace <- with(Data_alpha,aggregate(beta_trim,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(0,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(beta_trace_plus,na.rm = T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_minus,na.rm=T) + 
                                   colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(beta_trace_minus,na.rm=T) - 
                                                                                         colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_plus,na.rm=T) + 
                                   colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(beta_trace_plus,na.rm=T) - 
                                                                                       colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))

# Plot trace Beta experiment

plus_first <- unique(subset(Data_beta,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_beta,phase==0&condition=='minus')$sub)

alpha_trace <- with(Data_beta,aggregate(alpha_trim,by=list(sub=sub,trial=trial,condition=condition),mean))
alpha_trace_minus <- cast(subset(alpha_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
count_plus <- count_plus[2:length(count_plus)]
count_minus <- count_minus[2:length(count_minus)]

plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
     ylim=c(-10,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(alpha_trace_plus,na.rm=T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_minus,na.rm=T) + 
                                   colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(alpha_trace_minus,na.rm=T) - 
                                                                                          colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_plus,na.rm=T) + 
                                   colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(alpha_trace_plus,na.rm=T) - 
                                                                                        colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))


beta_trace <- with(Data_beta,aggregate(beta_trim,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(0,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(beta_trace_plus,na.rm = T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_minus,na.rm=T) + 
                                   colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(beta_trace_minus,na.rm=T) - 
                                                                                         colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_plus,na.rm=T) + 
                                   colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(beta_trace_plus,na.rm=T) - 
                                                                                       colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))

dev.off()



# Plot traces alpha learn -------------------------------------------------
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

jpeg(filename = "traces_alpha_learn.jpg",units = 'cm',width = 42,height = 30,res=300)
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

# Plot ExpA trace
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


# Plot ExpA predictions
plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)

conf_min <- with(subset(cj_pred_ma_alpha_learn,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_alpha_learn,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_alpha_learn,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_alpha_learn,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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

# Plot ExpB predictions 
plus_first <- unique(subset(Data_beta,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_beta,phase==0&condition=='minus')$sub)

conf_min <- with(subset(cj_pred_ma_alpha_learn,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_alpha_learn,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_alpha_learn,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_alpha_learn,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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

# Plot parameter traces 
par_higheta <- subset(par,eta>10000)

plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)

# Plot trace alpha experiment
alpha_trace <- with(Data_alpha,aggregate(alpha_alpha_learn,by=list(sub=sub,trial=trial,condition=condition),mean))
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
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_minus,na.rm=T) + 
                                   colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(alpha_trace_minus,na.rm=T) - 
                                                                                          colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_plus,na.rm=T) + 
                                   colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(alpha_trace_plus,na.rm=T) - 
                                                                                        colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))


beta_trace <- with(Data_alpha,aggregate(beta_alpha_learn,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(0,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(beta_trace_plus,na.rm = T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_minus,na.rm=T) + 
                                   colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(beta_trace_minus,na.rm=T) - 
                                                                                         colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_plus,na.rm=T) + 
                                   colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(beta_trace_plus,na.rm=T) - 
                                                                                       colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))

# Plot trace Beta experiment

plus_first <- unique(subset(Data_beta,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_beta,phase==0&condition=='minus')$sub)

alpha_trace <- with(Data_beta,aggregate(alpha_alpha_learn,by=list(sub=sub,trial=trial,condition=condition),mean))
alpha_trace_minus <- cast(subset(alpha_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
count_plus <- count_plus[2:length(count_plus)]
count_minus <- count_minus[2:length(count_minus)]

plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
     ylim=c(-10,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(alpha_trace_plus,na.rm=T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_minus,na.rm=T) + 
                                   colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(alpha_trace_minus,na.rm=T) - 
                                                                                          colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_plus,na.rm=T) + 
                                   colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(alpha_trace_plus,na.rm=T) - 
                                                                                        colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))


beta_trace <- with(Data_beta,aggregate(beta_alpha_learn,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(0,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(beta_trace_plus,na.rm = T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_minus,na.rm=T) + 
                                   colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(beta_trace_minus,na.rm=T) - 
                                                                                         colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_plus,na.rm=T) + 
                                   colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(beta_trace_plus,na.rm=T) - 
                                                                                       colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))

dev.off()



# Plot traces beta learn --------------------------------------------------
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

jpeg(filename = "traces_beta_learn.jpg",units = 'cm',width = 42,height = 30,res=300)
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

# Plot ExpA trace 
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


# Plot ExpA predictions
plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)

conf_min <- with(subset(cj_pred_ma_beta_learn,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_beta_learn,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_beta_learn,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_beta_learn,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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

# Plot ExpB predictions 
plus_first <- unique(subset(Data_beta,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_beta,phase==0&condition=='minus')$sub)

conf_min <- with(subset(cj_pred_ma_beta_learn,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_beta_learn,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_beta_learn,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_beta_learn,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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

# Plot parameter traces 
par_higheta <- subset(par,eta>10000)

plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)

# Plot trace alpha experiment
alpha_trace <- with(Data_alpha,aggregate(alpha_beta_learn,by=list(sub=sub,trial=trial,condition=condition),mean))
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
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_minus,na.rm=T) + 
                                   colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(alpha_trace_minus,na.rm=T) - 
                                                                                          colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_plus,na.rm=T) + 
                                   colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(alpha_trace_plus,na.rm=T) - 
                                                                                        colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))


beta_trace <- with(Data_alpha,aggregate(beta_beta_learn,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(0,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(beta_trace_plus,na.rm = T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_minus,na.rm=T) + 
                                   colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(beta_trace_minus,na.rm=T) - 
                                                                                         colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_plus,na.rm=T) + 
                                   colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(beta_trace_plus,na.rm=T) - 
                                                                                       colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))

# Plot trace Beta experiment

plus_first <- unique(subset(Data_beta,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_beta,phase==0&condition=='minus')$sub)

alpha_trace <- with(Data_beta,aggregate(alpha_beta_learn,by=list(sub=sub,trial=trial,condition=condition),mean))
alpha_trace_minus <- cast(subset(alpha_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
count_plus <- count_plus[2:length(count_plus)]
count_minus <- count_minus[2:length(count_minus)]

plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
     ylim=c(-10,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(alpha_trace_plus,na.rm=T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_minus,na.rm=T) + 
                                   colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(alpha_trace_minus,na.rm=T) - 
                                                                                          colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_plus,na.rm=T) + 
                                   colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(alpha_trace_plus,na.rm=T) - 
                                                                                        colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))


beta_trace <- with(Data_beta,aggregate(beta_beta_learn,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(0,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(beta_trace_plus,na.rm = T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_minus,na.rm=T) + 
                                   colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(beta_trace_minus,na.rm=T) - 
                                                                                         colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_plus,na.rm=T) + 
                                   colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(beta_trace_plus,na.rm=T) - 
                                                                                       colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))

dev.off()


# Plot traces both learn --------------------------------------------------
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

jpeg(filename = "traces_both_learn.jpg",units = 'cm',width = 42,height = 30,res=300)
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

# Plot ExpA trace 
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


# Plot ExpA predictions 
plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)

conf_min <- with(subset(cj_pred_ma_both_learn,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_both_learn,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_both_learn,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_both_learn,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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

# Plot ExpB predictions 
plus_first <- unique(subset(Data_beta,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_beta,phase==0&condition=='minus')$sub)

conf_min <- with(subset(cj_pred_ma_both_learn,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_both_learn,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_both_learn,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_both_learn,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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

# Plot parameter traces 
par_higheta <- subset(par,eta_a>10000)

plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)

# Plot trace alpha experiment
alpha_trace <- with(Data_alpha,aggregate(alpha_both_learn,by=list(sub=sub,trial=trial,condition=condition),mean))
alpha_trace_minus <- cast(subset(alpha_trace,sub %in% minus_first),sub ~trial, value = "x", fun.aggregate = mean)
alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first & !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
count_plus <- count_plus[2:length(count_plus)]
count_minus <- count_minus[2:length(count_minus)]

plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
     ylim = c(0,5),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(alpha_trace_plus,na.rm=T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_minus,na.rm=T) + 
                                   colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(alpha_trace_minus,na.rm=T) - 
                                                                                          colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_plus,na.rm=T) + 
                                   colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(alpha_trace_plus,na.rm=T) - 
                                                                                        colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))


beta_trace <- with(Data_alpha,aggregate(beta_both_learn,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(0,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(beta_trace_plus,na.rm = T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_minus,na.rm=T) + 
                                   colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(beta_trace_minus,na.rm=T) - 
                                                                                         colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_plus,na.rm=T) + 
                                   colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(beta_trace_plus,na.rm=T) - 
                                                                                       colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))

# Plot trace Beta experiment

plus_first <- unique(subset(Data_beta,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_beta,phase==0&condition=='minus')$sub)

alpha_trace <- with(Data_beta,aggregate(alpha_both_learn,by=list(sub=sub,trial=trial,condition=condition),mean))
alpha_trace_minus <- cast(subset(alpha_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
count_plus <- count_plus[2:length(count_plus)]
count_minus <- count_minus[2:length(count_minus)]

plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
     ylim=c(0,5),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(alpha_trace_plus,na.rm=T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_minus,na.rm=T) + 
                                   colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(alpha_trace_minus,na.rm=T) - 
                                                                                          colSds(alpha_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(alpha_trace_plus,na.rm=T) + 
                                   colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(alpha_trace_plus,na.rm=T) - 
                                                                                        colSds(alpha_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))


beta_trace <- with(Data_beta,aggregate(beta_both_learn,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first& !(sub %in% par_higheta$sub)),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(0,30),bty='n')
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
lines(colMeans(beta_trace_plus,na.rm = T),col=VERMILLION)
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_minus,na.rm=T) + 
                                   colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus),(colMeans(beta_trace_minus,na.rm=T) - 
                                                                                         colSds(beta_trace_minus,na.rm=T)/sqrt(count_minus))[Ntrials:1]),
        border=F,col=rgb(0,114,178,51,maxColorValue = 255))
polygon(c(1:Ntrials,Ntrials:1),c(colMeans(beta_trace_plus,na.rm=T) + 
                                   colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus),(colMeans(beta_trace_plus,na.rm=T) - 
                                                                                       colSds(beta_trace_plus,na.rm=T)/sqrt(count_plus))[Ntrials:1]),
        border=F,col=rgb(213,94,0,51,maxColorValue = 255))

dev.off()


# Trace per sub - NEED TO CHOOSE MODEL --------------------------------------------------------
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

# Model comparison plots --------------------------------------------------
# Plot parameters and cost distributions for all models
ggplot(par, aes(x = a0, colour = model, label = model)) +
  geom_density() +
  theme_bw()
ggplot(par, aes(x = b0, colour = model, label = model)) +
  geom_density() +
  theme_bw()
ggplot(par, aes(x = eta_a, colour = model, label = model)) +
  geom_density() +
  theme_bw()
ggplot(subset(par,model != "alpha"), aes(x = eta_b, colour = model, label = model)) +
  geom_density() +
  theme_bw()
ggplot(par, aes(x = cost_ldc, colour = model, label = model)) +
  geom_density() +
  theme_bw()

# Correlation btw both and beta model
par_beta <- subset(par,model=='beta')
par_both <- subset(par,model=='both')
scatterplot(par_beta$a0 ~ par_both$a0, bty = 'n',regLine=T,boxplots=F,smooth=F,
            xlab="Both learn",ylab = "Beta learn",
            main = paste("a0, r =",round(cor(par_beta$a0, par_both$a0),3)))
scatterplot(par_beta$b0 ~ par_both$b0, bty = 'n',regLine=T,boxplots=F,smooth=F,
            xlab="Both learn",ylab = "Beta learn",
            main = paste("b0, r =",round(cor(par_beta$b0, par_both$b0),3)))
scatterplot(par_beta$eta_b ~ par_both$eta_b, bty = 'n',regLine=T,boxplots=F,smooth=F,
            xlab="Both learn",ylab = "Beta learn",
            main = paste("learning rate, r =",round(cor(par_beta$eta_b, par_both$eta_b),3)))
cor.test(par_beta$a0, par_both$a0)
cor.test(par_beta$b0, par_both$b0)
cor.test(par_beta$eta_b, par_both$eta_b)

hist(par_both$eta_a,bty='n',main="2 Learning rates model",xlab="Alpha learning rate")