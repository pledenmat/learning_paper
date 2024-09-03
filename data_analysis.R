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
# library(hrbrthemes)
library(ggplot2)
library(geomtextpath)
sourceCpp("ldc_train.cpp")
source("ldc_nn_functions.R")
plots <- F
stat_tests <- F

# Function & plot parameters ----------------------------------------------------------------
max_count <- function(data){
  return(max(table(data)))
}

se <- function(x,na.rm=F) sd(x,na.rm=na.rm)/sqrt(length(x))

error.bar <- function(x, y, upper, lower=upper, length=0,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length,...)
}

se <- function(x,na.rm=F) sd(x,na.rm=na.rm)/sqrt(length(x))

#' Function to get the right font size in points
cex_size <- function(size,cex.layout) {
  return(size/(par()$ps*cex.layout))
}
### Adjust sizes and positions
cex.layout <- .66

title_line <- -2
cex.title <- 3
cex.lab <- 3
cex.axis <- 2
cex.legend <- 2
# Text
cex.legend <- cex_size(8,cex.layout)
cex.phase <- cex_size(11,cex.layout)*cex.layout
cex.main <- cex_size(12,cex.layout)
cex.axis <- cex_size(8,cex.layout) 
cex.lab <- cex_size(10,cex.layout)*cex.layout
cex.lab.rel <- cex_size(10,cex.layout) # In a big layout, cex for mtext and cex.lab differ
cex.trialtype <- cex_size(8,cex.layout)
# Lines and dots
lwd.dat <- 1.5
cex.datdot <- 1.125
cex.legend.square <- 3
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
                  c("sub","age","gender","manip","handedness","questionID","questionResp")]

# Remove questionnaire from main data frame
Data <- Data[complete.cases(Data$block),]
rm_col <- c("File","questionID","questionResp","X","?..sub")
Data <- Data[,-which(names(Data) %in% rm_col)]


subs <- unique(Data$sub); Nsub <- length(subs)
Ntrials <- dim(Data)[1]/Nsub

Data$response <- 0
Data[Data$resp=="['c']",'response'] <- 1

Data[Data$manip=='alfa','manip'] <- 'alpha'


# Demographics
age_alpha <- with(subset(Data,manip=='alpha'),aggregate(age,by=list(sub),mean))
summary(age_alpha$x)
gender_alpha <- table(subset(Data,manip=='alpha')$gender)/756

age_beta <- with(subset(Data,manip=='beta'),aggregate(age,by=list(sub),mean))
summary(age_beta$x)
sd(age_beta$x)
gender_beta <- table(subset(Data,manip=='beta')$gender)/756

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

table(unique(exclusion) %in% subset(Data,manip=='beta')$sub)

#' Filter out participants who reported only one confidence level 
#' more than 85% of the time
conf_count <- with(Data,aggregate(cj,by=list(sub=sub),max_count))
conf_count$x <- conf_count$x/Ntrials
exclusion <- c(exclusion, unique(conf_count[conf_count$x>.9,"sub"]))

table(unique(exclusion) %in% subset(Data,manip=='beta')$sub)

Data <- subset(Data,!(sub %in% exclusion))
questions <- subset(questions,!(sub %in% exclusion))

# Convert RTs to seconds
Data$rt <- Data$rt/1000
Data$RTconf <- Data$RTconf/1000

Data$response[Data$response==0] <- -1

Data$condition <- as.factor(Data$condition)
Data$difflevel <- as.factor(Data$difflevel)
# Data$sub <- as.factor(Data$sub)
subs <- unique(Data$sub); Nsub <- length(subs)

names(Data)[1] <- "trial"

# Confidence and feedback on scale from 0 to 1
Data$cj <- Data$cj/6
Data$fb <- Data$fb/100

Data$rt2 <- Data$rt + Data$RTconf # Total RT

# Separate participants according to which condition was first for plotting
Data$order <- "plus_first"
minus_first <- unique(subset(Data,phase==0&condition=='minus')$sub)
Data[Data$sub %in% minus_first,"order"] <- "minus_first"

setwd("..")
Nskip <- 0
Data <- subset(Data,RTconf<5 & rt<5 & rt>.2 & trial>=Nskip)
Ntrials <- Ntrials - Nskip # Because we now skip the first 5 trials
# write.csv(Data,file = "alternating_fb_mod_trim_skip.csv",row.names = F)
write.csv(Data,file = "alternating_fb_mod_trim.csv",row.names = F)
# write.csv(Data,file = "alternating_fb_mod.csv",row.names = F)

Nphase_trial <- length(unique(Data$withinphasetrial))
Nphase_block <- 18
Data$phase_block <-  Data$withinphasetrial %/% (Nphase_trial/Nphase_block)
Data$phase_block <- as.factor(Data$phase_block)

# Data$cor <- as.factor(Data$cor)

Data$cj_integer <- Data$cj*6

Data_alpha <- subset(Data,manip=='alpha')
Data_beta <- subset(Data,manip=='beta')
Nalpha <- length(unique(Data_alpha$sub))
Nbeta <- length(unique(Data_beta$sub))

# Behavior analysis - Dynamics -----------------------------------------------------------
if (stat_tests) {
  control <- lmerControl(optimizer = "bobyqa")
  glmercontrol <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5))
  
  # Set contrast coding
  options(contrasts=c("contr.sum","contr.poly"))
  
  ## Replication of previous experiments (static effects)
  # RT
  rt.int.alpha <- lmer(rt~condition*difflevel*withinphasetrial + (1|sub),data = subset(Data_alpha,cor==1),REML = F,control = control)
  rt.cond.alpha <- lmer(rt~condition*difflevel*withinphasetrial + (condition|sub),data = subset(Data_alpha,cor==1),REML = F,control = control)
  anova(rt.int.alpha,rt.cond.alpha)
  # Singular fit
  rt.cond.diff.alpha <- lmer(rt~condition*difflevel*withinphasetrial + (condition+difflevel|sub),data = subset(Data_alpha,cor==1),REML = F,control = control)
  anova(rt.cond.alpha)
  
  rt.int.beta <- lmer(rt~condition*difflevel*withinphasetrial + (1|sub),data = subset(Data_beta,cor==1),REML = F,control = control)
  rt.cond.beta <- lmer(rt~condition*difflevel*withinphasetrial + (condition|sub),data = subset(Data_beta,cor==1),REML = F,control = control)
  anova(rt.int.beta,rt.cond.beta)
  # Singular fit
  rt.cond.diff.beta <- lmer(rt~condition*difflevel*withinphasetrial + (condition+difflevel|sub),data = subset(Data_beta,cor==1),REML = F,control = control)
  anova(rt.cond.beta)
  
  # Accuracy
  acc.int.alpha <- glmer(cor~condition*difflevel*withinphasetrial + (1|sub),data = Data_alpha,family = binomial, control = glmercontrol)
  acc.cond.alpha <- glmer(cor~condition*difflevel*withinphasetrial + (condition|sub),data = Data_alpha,family = binomial, control = glmercontrol)
  anova(acc.int.alpha,acc.cond.alpha)
  acc.cond.diff.alpha <- glmer(cor~condition*difflevel*withinphasetrial + (condition+difflevel|sub),data = Data_alpha,family = binomial)
  Anova(acc.cond.alpha)
  
  acc.int.beta <- glmer(cor~condition*difflevel*withinphasetrial + (1|sub),data = Data_beta,family = binomial, control = glmercontrol)
  acc.cond.beta <- glmer(cor~condition*difflevel*withinphasetrial + (condition|sub),data = Data_beta,family = binomial, control = glmercontrol)
  anova(acc.int.beta,acc.cond.beta)
  acc.cond.diff.beta <- glmer(cor~condition*difflevel*withinphasetrial + (condition+difflevel|sub),data = Data_beta,family = binomial)
  Anova(acc.cond.beta)
  
  # Alpha
  cj.int.alpha <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (1|sub),data = Data_alpha,REML = F, ,control = control); 
  cj.cond.alpha <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (condition|sub),data = Data_alpha,REML = F, ,control = control); 
  anova(m.int,m.cond)
  cj.cond.acc.alpha <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (cor + condition|sub),data = Data_alpha, REML = F, ,control = control)
  anova(m.cond,m.cond.acc)
  # Singular fit
  # cj.cond.acc.difflevel.alpha <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (cor + condition + difflevel|sub),
  #                      data = Data_alpha, REML = F,control = control)
  cj.cond.acc.interaction.alpha <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (cor + condition + condition:cor|sub),
                                        data = Data_alpha, REML = F,control = control)
  plot(resid(cj.cond.acc.interaction.alpha),Data_alpha$cj) #Linearity
  leveneTest(residuals(cj.cond.acc.interaction.alpha) ~ Data_alpha$cor*Data_alpha$condition*Data_alpha$difflevel) #Homogeneity of variance
  qqmath(cj.cond.acc.interaction.alpha) #Normality
  anova(cj.cond.acc.interaction.alpha) #Results
  cj.cond.interaction.alpha.cor <- lmer(cj ~ condition*difflevel*withinphasetrial + (condition|sub),
                                        data = subset(Data_alpha,cor==1), REML = F,control = control)
  anova(cj.cond.interaction.alpha.cor)
  cj.cond.interaction.alpha.err <- lmer(cj ~ condition*difflevel*withinphasetrial + (condition|sub),
                                        data = subset(Data_alpha,cor==0), REML = F,control = control)
  anova(cj.cond.interaction.alpha.err)
  
  # Beta
  cj.int.beta <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (1|sub),data = Data_beta,REML = F, ,control = control); 
  cj.cond.beta <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (condition|sub),data = Data_beta,REML = F, ,control = control); 
  anova(cj.int.beta,cj.cond.beta)
  cj.cond.acc.beta <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (cor + condition|sub),data = Data_beta, REML = F, ,control = control)
  anova(cj.cond.beta,cj.cond.acc.beta)
  # Singular fit
  # cj.cond.acc.difflevel.beta <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (cor + condition + difflevel|sub),
  #                      data = Data_beta, REML = F,control = control)
  cj.cond.acc.interaction.beta <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (cor + condition + condition:cor|sub),
                                       data = Data_beta, REML = F,control = control)
  anova(cj.cond.acc.beta,cj.cond.acc.interaction.beta)
  plot(resid(cj.cond.acc.interaction.beta),Data_beta$cj) #Linearity
  leveneTest(residuals(cj.cond.acc.interaction.beta) ~ Data_beta$cor*Data_beta$condition*Data_beta$difflevel) #Homogeneity of variance
  qqmath(cj.cond.acc.interaction.beta) #Normality
  vif(cj.cond.acc.interaction.beta) # Multicollinearity
  anova(cj.cond.acc.interaction.beta) #Results  
}


# Behavior analysis - Static ----------------------------------------------
if (stat_tests) {
  ## Replication of previous experiments (static effects)
  # RT
  rt.int.alpha.static <- lmer(rt~condition*difflevel + (1|sub),data = subset(Data_alpha,cor==1),REML = F,control = control)
  rt.cond.alpha.static <- lmer(rt~condition*difflevel + (condition|sub),data = subset(Data_alpha,cor==1),REML = F,control = control)
  anova(rt.int.alpha.static,rt.cond.alpha.static)
  # Singular fit
  rt.cond.diff.alpha.static <- lmer(rt~condition*difflevel + (condition+difflevel|sub),data = subset(Data_alpha,cor==1),REML = F,control = control)
  anova(rt.cond.alpha.static)
  
  rt.int.beta.static <- lmer(rt~condition*difflevel + (1|sub),data = subset(Data_beta,cor==1),REML = F,control = control)
  rt.cond.beta.static <- lmer(rt~condition*difflevel + (condition|sub),data = subset(Data_beta,cor==1),REML = F,control = control)
  anova(rt.int.beta.static,rt.cond.beta.static)
  # Singular fit
  rt.cond.diff.beta.static <- lmer(rt~condition*difflevel + (condition+difflevel|sub),data = subset(Data_beta,cor==1),REML = F,control = control)
  anova(rt.cond.beta.static)
  
  # Accuracy
  acc.int.alpha.static <- glmer(cor~condition*difflevel + (1|sub),data = Data_alpha,family = binomial, control = glmercontrol)
  acc.cond.alpha.static <- glmer(cor~condition*difflevel + (condition|sub),data = Data_alpha,family = binomial, control = glmercontrol)
  anova(acc.int.alpha.static,acc.cond.alpha.static)
  acc.cond.diff.alpha.static <- glmer(cor~condition*difflevel + (condition+difflevel|sub),data = Data_alpha,family = binomial)
  Anova(acc.cond.alpha.static)
  
  acc.int.beta.static <- glmer(cor~condition*difflevel + (1|sub),data = Data_beta,family = binomial, control = glmercontrol)
  acc.cond.beta.static <- glmer(cor~condition*difflevel + (condition|sub),data = Data_beta,family = binomial, control = glmercontrol)
  anova(acc.int.beta.static,acc.cond.beta.static)
  acc.cond.diff.beta.static <- glmer(cor~condition*difflevel + (condition+difflevel|sub),data = Data_beta,family = binomial)
  Anova(acc.cond.beta.static)
  
  # Alpha
  cj.int.alpha.static <- lmer(cj ~ condition*cor*difflevel + (1|sub),data = Data_alpha,REML = F, ,control = control); 
  cj.cond.alpha.static <- lmer(cj ~ condition*cor*difflevel + (condition|sub),data = Data_alpha,REML = F, ,control = control); 
  anova(m.int,m.cond)
  cj.cond.acc.alpha.static <- lmer(cj ~ condition*cor*difflevel + (cor + condition|sub),data = Data_alpha, REML = F, ,control = control)
  anova(m.cond,m.cond.acc)
  # Singular fit
  # cj.cond.acc.difflevel.alpha.static <- lmer(cj ~ condition*cor*difflevel + (cor + condition + difflevel|sub),
  #                      data = Data_alpha, REML = F,control = control)
  cj.cond.acc.interaction.alpha.static <- lmer(cj ~ condition*cor*difflevel + (cor + condition + condition:cor|sub),
                                               data = Data_alpha, REML = F,control = control)
  plot(resid(cj.cond.acc.interaction.alpha.static),Data_alpha$cj) #Linearity
  leveneTest(residuals(cj.cond.acc.interaction.alpha.static) ~ Data_alpha$cor*Data_alpha$condition*Data_alpha$difflevel) #Homogeneity of variance
  qqmath(cj.cond.acc.interaction.alpha.static) #Normality
  vif(cj.cond.acc.interaction.alpha.static) # Multicollinearity
  anova(cj.cond.acc.interaction.alpha.static) #Results
  
  cj.cond.static.correct <- lmer(cj ~ condition*difflevel + (condition|sub),
                                 data = subset(Data_alpha,cor==1), REML = F,control = control)
  anova(cj.cond.static.correct)
  cj.cond.static.error <- lmer(cj ~ condition*difflevel + (condition|sub),
                               data = subset(Data_alpha,cor==0), REML = F,control = control)
  anova(cj.cond.static.error)
  
  # Beta
  cj.int.beta.static <- lmer(cj ~ condition*cor*difflevel + (1|sub),data = Data_beta,REML = F, ,control = control); 
  cj.cond.beta.static <- lmer(cj ~ condition*cor*difflevel + (condition|sub),data = Data_beta,REML = F, ,control = control); 
  anova(cj.int.beta.static,cj.cond.beta.static)
  cj.cond.acc.beta.static <- lmer(cj ~ condition*cor*difflevel + (cor + condition|sub),data = Data_beta, REML = F, ,control = control)
  anova(cj.cond.beta.static,cj.cond.acc.beta.static)
  # Singular fit
  # cj.cond.acc.difflevel.beta.static <- lmer(cj ~ condition*cor*difflevel + (cor + condition + difflevel|sub),
  #                      data = Data_beta, REML = F,control = control)
  cj.cond.acc.interaction.beta.static <- lmer(cj ~ condition*cor*difflevel + (cor + condition + condition:cor|sub),
                                              data = Data_beta, REML = F,control = control)
  anova(cj.cond.acc.beta.static,cj.cond.acc.interaction.beta.static)
  plot(resid(cj.cond.acc.interaction.beta.static),Data_beta$cj) #Linearity
  leveneTest(residuals(cj.cond.acc.interaction.beta.static) ~ Data_beta$cor*Data_beta$condition*Data_beta$difflevel) #Homogeneity of variance
  qqmath(cj.cond.acc.interaction.beta.static) #Normality
  anova(cj.cond.acc.interaction.beta.static) #Results  
  post_hoc <- emmeans(cj.cond.acc.interaction.beta.static, ~ difflevel)
  pairs(post_hoc)
}

# Plot Feedback presented ------------------------------------------------
width <- 16 # Plot size expressed in cm
height <- 10

se <- function(x,na.rm=F) sd(x,na.rm=na.rm)/sqrt(length(x))

title_line <- -2
cex.title <- 3
cex.lab <- 3
cex.lab.rel <- 3*.83
cex.axis <- 2
cex.legend <- 2
cex.datdot <- 2
lwd.dat <- 2
if (max(Data_alpha$fb)<=1) {
  Data_alpha$fb <- Data_alpha$fb*100
  Data_beta$fb <- Data_beta$fb*100
}
if (plots) {
  go_to("plots")
  go_to("paper")
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
  
  jpeg(filename = "feedback_presented.jpg",height = 16,width = 32,units = 'cm', res = 600)
  par(mfrow=c(1,2))
  stripchart(xminus_cor,ylim=c(0,100), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
             yaxt = 'n',xlab="",ylab = "")
  title(ylab = "Feedback", xlab = "Trial difficulty", line = 2.5,cex.lab=cex.lab.rel)
  axis(1,at=0:(n-1),labels=names(xminus_cor), cex.axis=cex.axis);
  axis(2, seq(0,100,20), cex.axis=cex.axis)
  means <- sapply(xminus_cor, mean)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=BLUE,lwd=lwd.dat,lty = "dashed")
  error.bar(0:(n-1),means,colSds(as.matrix(xminus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=BLUE)
  means <- sapply(xplus_cor, mean,na.rm=T)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=VERMILLION,lwd=lwd.dat,lty = "dashed")
  error.bar(0:(n-1),means,colSds(as.matrix(xplus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=VERMILLION)
  
  means <- sapply(xminus_err, mean, na.rm=T)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=BLUE,lwd=lwd.dat,lty = "dotted")
  error.bar(0:(n-1),means,colSds(as.matrix(xminus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=BLUE)
  means <- sapply(xplus_err, mean,na.rm=T)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=VERMILLION,lwd=lwd.dat,lty = "dotted")
  error.bar(0:(n-1),means,colSds(as.matrix(xplus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=VERMILLION)
  
  legend("bottomleft",legend=c("Correct trials","Error trials"),
         title = NULL,lty=c("dashed","dotted"),bty = "n",inset=0,
         cex = cex.legend, lwd=lwd.dat,seg.len=1.5)
  legend("bottomright",legend=c("Minus","Plus"),
         title = NULL,pch=rep(16,3),bty = "n",inset=0,
         cex = cex.legend,col=c(BLUE,VERMILLION))
  
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
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=BLUE,lwd=lwd.dat,lty = "dashed")
  error.bar(0:(n-1),means,colSds(as.matrix(xminus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=BLUE)
  means <- sapply(xplus_cor, mean,na.rm=T)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=VERMILLION,lwd=lwd.dat,lty = "dashed")
  error.bar(0:(n-1),means,colSds(as.matrix(xplus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=VERMILLION)
  
  means <- sapply(xminus_err, mean, na.rm=T)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=BLUE,lwd=lwd.dat,lty = "dotted")
  error.bar(0:(n-1),means,colSds(as.matrix(xminus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=BLUE)
  means <- sapply(xplus_err, mean,na.rm=T)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=VERMILLION,lwd=lwd.dat,lty = "dotted")
  error.bar(0:(n-1),means,colSds(as.matrix(xplus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=VERMILLION)
  
  legend("bottomleft",legend=c("Correct trials","Error trials"),
         title = NULL,lty=c("dashed","dotted"),bty = "n",inset=0,
         cex = cex.legend, lwd=lwd.dat,seg.len=1.5)
  dev.off()
  par(mfrow=c(1,1))
}
# Plot Accuracy and RT ----------------------------------------------------
if (plots) {
  if (is.factor(Data_alpha$cor)) {
    Data_alpha$cor <- as.numeric(Data_alpha$cor)-1
    Data_beta$cor <- as.numeric(Data_beta$cor)-1
  }
  go_to("plots")
  go_to("paper")
  ###
  # Experiment A
  ###
  N_temp <- length(unique(Data_alpha$sub))
  
  rtlow <- with(subset(Data_alpha,condition=="minus"),aggregate(rt,by=list(sub,difflevel),mean));
  names(rtlow) <- c('sub','difflevel','rt')
  rtlow <- cast(rtlow,sub~difflevel)
  
  rthigh <- with(subset(Data_alpha,condition=="plus"),aggregate(rt,by=list(sub,difflevel),mean));
  names(rthigh) <- c('sub','difflevel','rt')
  rthigh <- cast(rthigh,sub~difflevel)
  
  # Drop subject column
  rtlow <- rtlow[,c(2:4)]
  rthigh <- rthigh[,c(2:4)]
  n <- length(rtlow)
  
  rtlow <- rtlow[,c("hard","medium","easy")];
  rthigh <- rthigh[,c("hard","medium","easy")]
  
  corlow <- with(subset(Data_alpha,condition=="minus"),aggregate(cor,by=list(sub,difflevel),mean));
  names(corlow) <- c('sub','difflevel','cor')
  corlow <- cast(corlow,sub~difflevel)
  
  corhigh <- with(subset(Data_alpha,condition=="plus"),aggregate(cor,by=list(sub,difflevel),mean));
  names(corhigh) <- c('sub','difflevel','cor')
  corhigh <- cast(corhigh,sub~difflevel)
  
  # Drop subject column
  corlow <- corlow[,c(2:4)]
  corhigh <- corhigh[,c(2:4)]
  n <- length(corlow)
  
  corlow <- corlow[,c("hard","medium","easy")];
  corhigh <- corhigh[,c("hard","medium","easy")]
  
  jpeg(filename = "objective_performance.jpg",height = 32,width = 32,units = 'cm', res = 600)
  layout(matrix(c(1,2,3,4),ncol=2,byrow = F))
  stripchart(rtlow,ylim=c(0.8,1.1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
             yaxt = 'n',xlab="",ylab = "", main = "Alpha-manipulated feedback", cex.main = cex.title)
  title(ylab = "Reaction time (s)", xlab = "Trial difficulty", line = 2.5,cex.lab=cex.lab.rel)
  axis(1,at=0:(n-1),labels=names(rtlow), cex.axis=cex.axis);
  axis(2, seq(0.8,1.1,.1), cex.axis=cex.axis)
  means <- sapply(rtlow, mean)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=BLUE,lwd=lwd.dat)
  error.bar(0:(n-1),means,colSds(as.matrix(rtlow),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=BLUE)
  
  means <- sapply(rthigh, mean, na.rm=T)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=VERMILLION,lwd=lwd.dat)
  error.bar(0:(n-1),means,colSds(as.matrix(rthigh),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=VERMILLION)
  
  legend("top",legend=c("low","high"),horiz = T,
         title = NULL,pch=rep(16,3),bty = "n",inset=0,
         cex = cex.legend,col=c(BLUE,VERMILLION))
  
  stripchart(corlow,ylim=c(0.5,1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
             yaxt = 'n',xlab="",ylab = "")
  title(ylab = "Accuracy", xlab = "Trial difficulty", line = 2.5,cex.lab=cex.lab.rel)
  axis(1,at=0:(n-1),labels=names(corlow), cex.axis=cex.axis);
  axis(2, seq(0.5,1,.1), cex.axis=cex.axis)
  means <- sapply(corlow, mean)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=BLUE,lwd=lwd.dat)
  error.bar(0:(n-1),means,colSds(as.matrix(corlow),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=BLUE)
  
  means <- sapply(corhigh, mean, na.rm=T)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=VERMILLION,lwd=lwd.dat)
  error.bar(0:(n-1),means,colSds(as.matrix(corhigh),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=VERMILLION)
  
  # legend("bottomright",legend=c("low","high"),
  #        title = NULL,pch=rep(16,3),bty = "n",inset=0,
  #        cex = cex.legend,col=c(BLUE,VERMILLION))
  
  ###
  # Experiment B
  ###
  N_temp <- length(unique(Data_beta$sub))
  
  rtlow <- with(subset(Data_beta,condition=="minus"),aggregate(rt,by=list(sub,difflevel),mean));
  names(rtlow) <- c('sub','difflevel','rt')
  rtlow <- cast(rtlow,sub~difflevel)
  
  rthigh <- with(subset(Data_beta,condition=="plus"),aggregate(rt,by=list(sub,difflevel),mean));
  names(rthigh) <- c('sub','difflevel','rt')
  rthigh <- cast(rthigh,sub~difflevel)
  
  
  # Drop subject column
  rtlow <- rtlow[,c(2:4)]
  rthigh <- rthigh[,c(2:4)]
  n <- length(rthigh)
  
  rtlow <- rtlow[,c("hard","medium","easy")];
  rthigh <- rthigh[,c("hard","medium","easy")];
  
  corlow <- with(subset(Data_beta,condition=="minus"),aggregate(cor,by=list(sub,difflevel),mean));
  names(corlow) <- c('sub','difflevel','cor')
  corlow <- cast(corlow,sub~difflevel)
  
  corhigh <- with(subset(Data_beta,condition=="plus"),aggregate(cor,by=list(sub,difflevel),mean));
  names(corhigh) <- c('sub','difflevel','cor')
  corhigh <- cast(corhigh,sub~difflevel)
  
  # Drop subject column
  corlow <- corlow[,c(2:4)]
  corhigh <- corhigh[,c(2:4)]
  n <- length(corlow)
  
  corlow <- corlow[,c("hard","medium","easy")];
  corhigh <- corhigh[,c("hard","medium","easy")]
  
  stripchart(rtlow,ylim=c(0.8,1.1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
             yaxt = 'n',xlab="",ylab = "",  main = "Beta-manipulated feedback", cex.main = cex.title)
  title(ylab = "Reaction time (s)", xlab = "Trial difficulty", line = 2.5,cex.lab=cex.lab.rel)
  axis(1,at=0:(n-1),labels=names(rtlow), cex.axis=cex.axis);
  axis(2, seq(0.8,1.1,.1), cex.axis=cex.axis)
  means <- sapply(rtlow, mean)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=BLUE,lwd=lwd.dat)
  error.bar(0:(n-1),means,colSds(as.matrix(rtlow),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=BLUE)
  
  means <- sapply(rthigh, mean, na.rm=T)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=VERMILLION,lwd=lwd.dat)
  error.bar(0:(n-1),means,colSds(as.matrix(rthigh),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=VERMILLION)
  
  legend("top",legend=c("low","high"),horiz=T,
         title = NULL,pch=rep(16,3),bty = "n",inset=0,
         cex = cex.legend,col=c(BLUE,VERMILLION))
  
  stripchart(corlow,ylim=c(0.5,1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
             yaxt = 'n',xlab="",ylab = "")
  title(ylab = "Accuracy", xlab = "Trial difficulty", line = 2.5,cex.lab=cex.lab.rel)
  axis(1,at=0:(n-1),labels=names(corlow), cex.axis=cex.axis);
  axis(2, seq(0.5,1,.1), cex.axis=cex.axis)
  means <- sapply(corlow, mean)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=BLUE,lwd=lwd.dat)
  error.bar(0:(n-1),means,colSds(as.matrix(corlow),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=BLUE)
  
  means <- sapply(corhigh, mean, na.rm=T)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=VERMILLION,lwd=lwd.dat)
  error.bar(0:(n-1),means,colSds(as.matrix(corhigh),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=VERMILLION)
  
  dev.off()
  par(mfrow=c(1,1))
}




# Plot static confidence - WIP --------------------------------------------

if (plots) {
  if (is.factor(Data_alpha$cor)) {
    Data_alpha$cor <- as.numeric(Data_alpha$cor)-1
    Data_beta$cor <- as.numeric(Data_beta$cor)-1
  }
  go_to("plots")
  go_to("paper")
  ###
  # Experiment A
  ###
  N_temp <- length(unique(Data_alpha$sub))
  
  cjlow <- with(subset(Data_alpha,condition=="minus"),aggregate(cj,by=list(sub,difflevel,cor),mean));
  names(cjlow) <- c('sub','difflevel','cor','rt')
  cjlow_err <- subset(cjlow,cor==0)
  cjlow_cor <- subset(cjlow,cor==1)
  cjlow_err <- cast(cjlow_err,sub~difflevel)
  cjlow_cor <- cast(cjlow_cor,sub~difflevel)
  
  cjhigh <- with(subset(Data_alpha,condition=="plus"),aggregate(cj,by=list(sub,difflevel,cor),mean));
  names(cjhigh) <- c('sub','difflevel','cor','rt')
  cjhigh_err <- subset(cjhigh,cor==0)
  cjhigh_cor <- subset(cjhigh,cor==1)
  cjhigh_err <- cast(cjhigh_err,sub~difflevel)
  cjhigh_cor <- cast(cjhigh_cor,sub~difflevel)
  
  # Drop subject column
  cjlow_err <- cjlow_err[,c(2:4)]
  cjlow_cor <- cjlow_cor[,c(2:4)]
  cjhigh_err <- cjhigh_err[,c(2:4)]
  cjhigh_cor <- cjhigh_cor[,c(2:4)]
  n <- length(cjlow_err)
  
  cjlow_err <- cjlow_err[,c("hard","medium","easy")];
  cjlow_cor <- cjlow_cor[,c("hard","medium","easy")]
  cjhigh_err <- cjhigh_err[,c("hard","medium","easy")]
  cjhigh_cor <- cjhigh_cor[,c("hard","medium","easy")]
  
  
  jpeg(filename = "static_confidence.jpg",height = 16,width = 32,units = 'cm', res = 600)
  par(mfrow=c(1,2))
  stripchart(cjlow_err,ylim=c(.6,1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
             yaxt = 'n',xlab="",ylab = "", main = "Alpha-manipulated feedback", cex.main = cex.title)
  title(ylab = "Confidence", xlab = "Trial difficulty", line = 2.5,cex.lab=cex.lab.rel)
  axis(1,at=0:(n-1),labels=names(cjlow_err), cex.axis=cex.axis);
  axis(2, seq(.6,1,.5), cex.axis=cex.axis)
  means <- sapply(cjlow_err, mean)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=BLUE,lwd=lwd.dat,lty='dotted')
  error.bar(0:(n-1),means,colSds(as.matrix(cjlow_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=BLUE,lty='dotted')
  
  means <- sapply(cjhigh_err, mean, na.rm=T)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=VERMILLION,lwd=lwd.dat,lty='dotted')
  error.bar(0:(n-1),means,colSds(as.matrix(cjhigh_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=VERMILLION,lty='dotted')
  
  means <- sapply(cjlow_cor, mean)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=BLUE,lwd=lwd.dat,lty='dashed')
  error.bar(0:(n-1),means,colSds(as.matrix(cjlow_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=BLUE,lty='dashed')
  
  means <- sapply(cjhigh_cor, mean, na.rm=T)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=VERMILLION,lwd=lwd.dat,lty='dashed')
  error.bar(0:(n-1),means,colSds(as.matrix(cjhigh_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=VERMILLION,lty='dashed')
  
  legend("top",legend=c("low","high"),horiz = T,
         title = NULL,pch=rep(16,3),bty = "n",inset=0,
         cex = cex.legend,col=c(BLUE,VERMILLION))
  legend("bottomleft",legend=c("Correct trials","Error trials"),
         title = NULL,lty=c("dashed","dotted"),bty = "n",inset=0,
         cex = cex.legend, lwd=lwd.dat,seg.len=1.5)
  
  ###
  # Experiment B
  ###
  N_temp <- length(unique(Data_beta$sub))
  
  cjlow <- with(subset(Data_beta,condition=="minus"),aggregate(cj,by=list(sub,difflevel,cor),mean));
  names(cjlow) <- c('sub','difflevel','cor','rt')
  cjlow_err <- subset(cjlow,cor==0)
  cjlow_cor <- subset(cjlow,cor==1)
  cjlow_err <- cast(cjlow_err,sub~difflevel)
  cjlow_cor <- cast(cjlow_cor,sub~difflevel)
  
  cjhigh <- with(subset(Data_beta,condition=="plus"),aggregate(cj,by=list(sub,difflevel,cor),mean));
  names(cjhigh) <- c('sub','difflevel','cor','rt')
  cjhigh_err <- subset(cjhigh,cor==0)
  cjhigh_cor <- subset(cjhigh,cor==1)
  cjhigh_err <- cast(cjhigh_err,sub~difflevel)
  cjhigh_cor <- cast(cjhigh_cor,sub~difflevel)
  
  # Drop subject column
  cjlow_err <- cjlow_err[,c(2:4)]
  cjlow_cor <- cjlow_cor[,c(2:4)]
  cjhigh_err <- cjhigh_err[,c(2:4)]
  cjhigh_cor <- cjhigh_cor[,c(2:4)]
  n <- length(cjlow_err)
  
  cjlow_err <- cjlow_err[,c("hard","medium","easy")];
  cjlow_cor <- cjlow_cor[,c("hard","medium","easy")]
  cjhigh_err <- cjhigh_err[,c("hard","medium","easy")]
  cjhigh_cor <- cjhigh_cor[,c("hard","medium","easy")]
  
  stripchart(cjlow_err,ylim=c(.6,1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
             yaxt = 'n',xlab="",ylab = "", main = "Beta-manipulated feedback", cex.main = cex.title)
  title(ylab = "Confidence", xlab = "Trial difficulty", line = 2.5,cex.lab=cex.lab.rel)
  axis(1,at=0:(n-1),labels=names(cjlow_err), cex.axis=cex.axis);
  axis(2, seq(.6,1,.5), cex.axis=cex.axis)
  means <- sapply(cjlow_err, mean)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=BLUE,lwd=lwd.dat,lty='dotted')
  error.bar(0:(n-1),means,colSds(as.matrix(cjlow_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=BLUE,lty='dotted')
  
  means <- sapply(cjhigh_err, mean, na.rm=T)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=VERMILLION,lwd=lwd.dat,lty='dotted')
  error.bar(0:(n-1),means,colSds(as.matrix(cjhigh_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=VERMILLION,lty='dotted')
  
  means <- sapply(cjlow_cor, mean)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=BLUE,lwd=lwd.dat,lty='dashed')
  error.bar(0:(n-1),means,colSds(as.matrix(cjlow_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=BLUE,lty='dashed')
  
  means <- sapply(cjhigh_cor, mean, na.rm=T)
  lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=VERMILLION,lwd=lwd.dat,lty='dashed')
  error.bar(0:(n-1),means,colSds(as.matrix(cjhigh_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=VERMILLION,lty='dashed')
  
  legend("top",legend=c("low","high"),horiz = T,
         title = NULL,pch=rep(16,3),bty = "n",inset=0,
         cex = cex.legend,col=c(BLUE,VERMILLION))
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
models <- c("no","alpha","beta","both")

totlen <- length(difflevels)*length(subs)*length(models)

# DDM parameters + hyperparameters + error
par <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                  bound = NA, drift = NA, ter = NA,cost_ddm = NA, 
                  sub = rep(subs, each = length(difflevels)*length(models)),
                  difflevel = rep(difflevels, length.out = totlen),
                  manip = NA, model = rep(models, each = length(difflevels)) )

Data$drift <- NA

go_to("fit")
go_to("alternating_fb")
if (!file.exists("anal_sim.Rdata")) {
  anal_sim <- Data[rep(seq_len(nrow(Data)), each=Nsim*length(models)), c('trial','withinphasetrial','sub','condition','cor')]
  anal_sim$alpha <- NA
  anal_sim$beta <- NA
  anal_sim$sim <- 1:Nsim
  anal_sim$model <- rep(models,each=Nsim)
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
  
  par[par$sub==subs[s],"bound"] <- ddm.results$optim$bestmem[1]
  par[par$sub==subs[s],"ter"] <- ddm.results$optim$bestmem[2]
  par[par$sub==subs[s],"z"] <- ddm.results$optim$bestmem[3]
  par[par$sub==subs[s],"vratio"] <- ddm.results$optim$bestmem[4]
  par[par$sub==subs[s]&par$difflevel==difflevels[1],"drift"] <- ddm.results$optim$bestmem[5]
  par[par$sub==subs[s]&par$difflevel==difflevels[2],"drift"] <- ddm.results$optim$bestmem[6]
  par[par$sub==subs[s]&par$difflevel==difflevels[3],"drift"] <- ddm.results$optim$bestmem[7]
  par[par$sub==subs[s],"cost_ddm"] <- ddm.results$optim$bestval
  
  #Add estimated drift to dataset
  for (diff in difflevels) {
    Data$drift[Data$difflevel==diff&Data$sub==subs[s]] <- unique(par[par$sub==subs[s]&par$difflevel==diff,"drift"])
  }
  temp_dat <- subset(Data,sub==subs[s])
  
  for (model in models) {
    ldc_file <- paste0('ldc_nn/trim/',model,'_learn/ldcfit_',subs[s],'.Rdata')
    if (file.exists(ldc_file)) {
      load(ldc_file)
      par[par$sub==subs[s]&par$model==model,"a0"] <- ldc.results$optim$bestmem[1]
      par[par$sub==subs[s]&par$model==model,"b0"] <- ldc.results$optim$bestmem[2]
      par[par$sub==subs[s]&par$model==model,"eta_a"] <- ldc.results$optim$bestmem[4]
      par[par$sub==subs[s]&par$model==model,"eta_b"] <- ldc.results$optim$bestmem[5]
      par[par$sub==subs[s]&par$model==model,"cost_ldc"] <- ldc.results$optim$bestval
      
      #' Generate model predictions multiple times to account for the stochastic
      #' nature of estimating single trial accumulated evidence 
      if (!file.exists("anal_sim2.Rdata")) {
        for (i in 1:Nsim) {
          results <-
            ldc.nn.fit.w(params=c(mean(par[par$sub==subs[s]&par$model==model,"a0"]),
                                  mean(par[par$sub==subs[s]&par$model==model,"b0"]),1,
                                  mean(par[par$sub==subs[s]&par$model==model,"eta_a"]),
                                  mean(par[par$sub==subs[s]&par$model==model,"eta_b"])),
                         ddm_params = ddm_params,
                         obs=temp_dat,returnFit = F,eta_sep=T,
                         Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                         dt = dt, sigma = sigma)
          anal_sim[anal_sim$sim==i&anal_sim$sub==subs[s]&anal_sim$model==model ,'cj'] <- results$pred
          anal_sim[anal_sim$sim==i&anal_sim$sub==subs[s]&anal_sim$model==model ,'alpha'] <- results$trace[,1]
          anal_sim[anal_sim$sim==i&anal_sim$sub==subs[s]&anal_sim$model==model ,'beta'] <- results$trace[,2]  
        }
      }
    }
  }
  
  par[par$sub==subs[s],"manip"] <- unique(temp_dat$manip)
}
# Save model simulations
if (!file.exists("anal_sim.Rdata")) {
  save(anal_sim,file="anal_sim_sim.Rdata")
} else {
  load("anal_sim.Rdata")
}

#' Add confidence prediction and parameter trace from each model to the empirical data frame
for (m in models) {
  print(m)
  cj_pred <- with(subset(anal_sim,model==m),aggregate(cj,by=list(trial,sub),mean))
  alpha <- with(subset(anal_sim,model==m),aggregate(alpha,by=list(trial,sub),mean))
  beta <- with(subset(anal_sim,model==m),aggregate(beta,by=list(trial,sub),mean))
  names(cj_pred) <- c("trial","sub",paste0("cj_pred_",m,"_learn"))
  names(alpha) <- c("trial","sub",paste0("alpha_",m,"_learn"))
  names(beta) <- c("trial","sub",paste0("beta_",m,"_learn"))
  anal_sim_mean <- merge(merge(cj_pred,alpha),beta)
  Data <- merge(Data,anal_sim_mean)
  
}

Data_alpha <- subset(Data,manip=='alpha')
Data_beta <- subset(Data,manip=='beta')
# Model comparison --------------------------------------------------------
par$Ndata_point <-  round(nrow(Data)/Nsub)
par$Npar <- 3
par[par$model=="no",'Npar'] <- 2
par[par$model=="both",'Npar'] <- 4

bic_custom <- function(Residuals,k,n){
  return(log(n)*k+n*log(Residuals/n))
}

par$bic <- bic_custom(par$cost_ldc,par$Npar,par$Ndata_point)
mean_bic <- with(par,aggregate(bic,by=list(model=model,manip=manip),mean))
mean_bic$delta <- -99
mean_bic[mean_bic$manip=="alpha","delta"] <- 
  mean_bic[mean_bic$manip=="alpha",]$x - 
  min(mean_bic[mean_bic$manip=="alpha",]$x)
mean_bic[mean_bic$manip=="beta","delta"] <- 
  mean_bic[mean_bic$manip=="beta",]$x -
  min(mean_bic[mean_bic$manip=="beta",]$x)

# Plot differences in BIC
jpeg("model_comparison.jpg",width=10,height=10,units = 'cm',res=600)
bp <- barplot(subset(mean_bic,manip=='beta')$delta,ylab = expression(paste(Delta,"BIC")), xlab = "")
# Add model names under each bar
bp.labels <- subset(mean_bic,manip=='beta')$model
bp.labels[bp.labels=="alpha"] <- expression(paste(alpha,"-learning"))
bp.labels[bp.labels=="beta"] <- expression(paste(beta,"-learning"))
bp.labels[bp.labels=="both"] <- "Full learning"
bp.labels[bp.labels=="no"] <- "No learning"
text(x = bp, y = -0.5, labels = bp.labels, srt = 45, adj = 1, xpd = TRUE)
dev.off()
### Best model per participant
model_bic <- with(par,aggregate(bic,list(sub=sub,manip=manip,model=model),mean))
model_bic <- cast(model_bic, sub + manip ~ model)
models_ord <- names(model_bic)[3:6]
model_bic$best <- models_ord[apply(model_bic[,3:6],1,which.min)]
with(model_bic,aggregate(best,list(manip),table))
table(subset(model_bic,manip=='alpha')$best)
table(subset(model_bic,manip=='beta')$best)

sub_nolearn_best <- subset(model_bic,best=='no')$sub
sub_learn_best <- subset(model_bic,best!='no')$sub

model_bic_learn_vs_nolearn <- subset(model_bic,manip=='beta')[,c('sub','manip','both','no','best')]
model_bic_learn_vs_nolearn$best <- c("both","no")[apply(model_bic_learn_vs_nolearn[,3:4],1,which.min)]
table(model_bic_learn_vs_nolearn$best)

# Add the best model to the empirical data frame
if ("best" %in% names(Data)) {
  Data <- Data[,!(names(Data) == "best")]
}
Data <- merge(Data,model_bic[,c('sub','manip','best')])
Data$best <- as.factor(Data$best)

# Add a column with the predicted confidence from the best model per participant to the empirical data frame
Data$cj_pred_best <- NA
for (s in 1:Nsub) {
  tempdat <- subset(Data,sub==subs[s])
  if (unique(tempdat$best=="no")) {
    Data[Data$sub==subs[s]&Data$manip==tempdat$manip[1],'cj_pred_best'] <- tempdat$cj_pred_no_learn
  } else if (unique(tempdat$best=="alpha")) {
    Data[Data$sub==subs[s]&Data$manip==tempdat$manip[1],'cj_pred_best'] <- tempdat$cj_pred_alpha_learn
  } else if (unique(tempdat$best=="beta")) {
    Data[Data$sub==subs[s]&Data$manip==tempdat$manip[1],'cj_pred_best'] <- tempdat$cj_pred_beta_learn
  } else if (unique(tempdat$best=="both")) {
    Data[Data$sub==subs[s]&Data$manip==tempdat$manip[1],'cj_pred_best'] <- tempdat$cj_pred_both_learn
  }
}

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

pred_conf_sub_no_learn <- with(Data,aggregate(cj_pred_no_learn,by=list(trial,cor,sub),mean))
names(pred_conf_sub_no_learn) <- c("trial","cor","sub","cj")

trials <- data.frame(trial=rep((0:(Ntrials-1))+Nskip,each=2),
                     cor=c(0,1),sub=rep(subs,each=Ntrials*2))

cj_ma <- merge(trial_conf_sub,trials,all=T)
cj_pred_ma_beta_learn <- merge(pred_conf_sub_beta_learn,trials,all=T)
cj_pred_ma_alpha_learn <- merge(pred_conf_sub_alpha_learn,trials,all=T)
cj_pred_ma_both_learn <- merge(pred_conf_sub_both_learn,trials,all=T)
cj_pred_ma_no_learn <- merge(pred_conf_sub_no_learn,trials,all=T)


ma <- function(x,n,names){
  return(rollapply(x[,names], width=n, FUN=function(x) mean(x, na.rm=TRUE),partial=TRUE, align="center"))
}
for (s in subs) {
  print(s)
  cj_ma[cj_ma$sub==s&cj_ma$cor==0,"cj"] <- ma(subset(cj_ma,sub==s&cor==0),n_err,"cj")
  cj_ma[cj_ma$sub==s&cj_ma$cor==1,"cj"] <- ma(subset(cj_ma,sub==s&cor==1),n,"cj")
  cj_pred_ma_beta_learn[cj_pred_ma_beta_learn$sub==s&cj_pred_ma_beta_learn$cor==0,"cj"] <- ma(subset(cj_pred_ma_beta_learn,sub==s&cor==0),n_err,"cj")
  cj_pred_ma_beta_learn[cj_pred_ma_beta_learn$sub==s&cj_pred_ma_beta_learn$cor==1,"cj"] <- ma(subset(cj_pred_ma_beta_learn,sub==s&cor==1),n,"cj")
  cj_pred_ma_alpha_learn[cj_pred_ma_alpha_learn$sub==s&cj_pred_ma_alpha_learn$cor==0,"cj"] <- ma(subset(cj_pred_ma_alpha_learn,sub==s&cor==0),n_err,"cj")
  cj_pred_ma_alpha_learn[cj_pred_ma_alpha_learn$sub==s&cj_pred_ma_alpha_learn$cor==1,"cj"] <- ma(subset(cj_pred_ma_alpha_learn,sub==s&cor==1),n,"cj")
  cj_pred_ma_both_learn[cj_pred_ma_both_learn$sub==s&cj_pred_ma_both_learn$cor==0,"cj"] <- ma(subset(cj_pred_ma_both_learn,sub==s&cor==0),n_err,"cj")
  cj_pred_ma_both_learn[cj_pred_ma_both_learn$sub==s&cj_pred_ma_both_learn$cor==1,"cj"] <- ma(subset(cj_pred_ma_both_learn,sub==s&cor==1),n,"cj")
  cj_pred_ma_no_learn[cj_pred_ma_no_learn$sub==s&cj_pred_ma_no_learn$cor==0,"cj"] <- ma(subset(cj_pred_ma_no_learn,sub==s&cor==0),n_err,"cj")
  cj_pred_ma_no_learn[cj_pred_ma_no_learn$sub==s&cj_pred_ma_no_learn$cor==1,"cj"] <- ma(subset(cj_pred_ma_no_learn,sub==s&cor==1),n,"cj")
}

cj_pred_ma_beta_learn$model <- "beta"
cj_pred_ma_alpha_learn$model <- "alpha"
cj_pred_ma_no_learn$model <- "no"
cj_pred_ma_both_learn$model <- "both"
cj_pred <- rbind(cj_pred_ma_no_learn,cj_pred_ma_alpha_learn,
                 cj_pred_ma_beta_learn,cj_pred_ma_both_learn)
# Plot traces LR ---------------------------------------------------

width <- 16 # Plot size expressed in cm
height <- 10


go_to("plots")
go_to("alternating_fb")
go_to("trim")
for (m in models) {
  jpeg(filename = paste0("traces_",m,"_learn.jpg"),units = 'cm',width = 42,height = 30,res=300)
  # layout(matrix(c(1,1,3,3,2,2,4,4,5,6,7,8),ncol=3))
  layout(matrix(c(1,2,9,9,5,6,11,11,1,3,10,10,5,7,12,12,1,4,13,14,5,8,15,16),ncol=3),heights = c(.05,.05,.2,.2,.05,.05,.2,.2))
  par(mar=c(0,0,0,0))
  plot.new()
  title(cex.main=cex.title,line=title_line,main = expression(paste(alpha,"-Manipulated Feedback, empirical fit")))
  plot.new()
  title(cex.main=cex.title,line=title_line,main= expression("Empirical Data"))
  plot.new()
  title(cex.main=cex.title,line=title_line,main= expression("Model Fits"))
  plot.new()
  title(cex.main=cex.title,line=title_line,main= expression("Weight traces"))
  plot.new()
  title(cex.main=cex.title,line=title_line,main = expression(paste(beta,"-Manipulated Feedback, empirical fit")))
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
  
  conf_min <- with(subset(cj_pred,model==m&sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
  names(conf_min) <- c("trial","cor","cj")
  conf_min_se <- with(subset(cj_pred,model==m&sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
  names(conf_min_se) <- c("trial","cor","cj")
  conf_plus <- with(subset(cj_pred,model==m&sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
  names(conf_plus) <- c("trial","cor","cj")
  conf_plus_se <- with(subset(cj_pred,model==m&sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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
  
  conf_min <- with(subset(cj_pred,model==m&sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
  names(conf_min) <- c("trial","cor","cj")
  conf_min_se <- with(subset(cj_pred,model==m&sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
  names(conf_min_se) <- c("trial","cor","cj")
  conf_plus <- with(subset(cj_pred,model==m&sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
  names(conf_plus) <- c("trial","cor","cj")
  conf_plus_se <- with(subset(cj_pred,model==m&sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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
  plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
  minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)
  
  # Plot trace alpha experiment
  alpha_trace <- Data_alpha[,c(paste0("alpha_",m,"_learn"),"sub","trial","condition")]
  names(alpha_trace) <- c("alpha","sub","trial","condition")
  alpha_trace_minus <- cast(subset(alpha_trace,sub %in% minus_first),sub ~trial, value = "alpha", fun.aggregate = mean)
  alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first),sub~trial, value = "alpha", fun.aggregate = mean)
  count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
  count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
  count_plus <- count_plus[2:length(count_plus)]
  count_minus <- count_minus[2:length(count_minus)]
  
  plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
       ylim = alpha_range,bty='n')
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
  
  
  beta_trace <- Data_alpha[,c(paste0("beta_",m,"_learn"),"sub","trial","condition")]
  names(beta_trace) <- c("beta","sub","trial","condition")
  beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "beta", fun.aggregate = mean)
  beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first),sub~trial, value = "beta", fun.aggregate = mean)
  plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
       ylim=beta_range,bty='n')
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
  
  alpha_trace <- Data_beta[,c(paste0("alpha_",m,"_learn"),"sub","trial","condition")]
  names(alpha_trace) <- c("alpha","sub","trial","condition")
  alpha_trace_minus <- cast(subset(alpha_trace,sub %in% minus_first),sub~trial, value = "alpha", fun.aggregate = mean)
  alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first),sub~trial, value = "alpha", fun.aggregate = mean)
  count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
  count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
  count_plus <- count_plus[2:length(count_plus)]
  count_minus <- count_minus[2:length(count_minus)]
  
  plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
       ylim=alpha_range,bty='n')
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
  
  
  beta_trace <- Data_beta[,c(paste0("beta_",m,"_learn"),"sub","trial","condition")]
  names(beta_trace) <- c("beta","sub","trial","condition")
  beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "beta", fun.aggregate = mean)
  beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first),sub~trial, value = "beta", fun.aggregate = mean)
  plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
       ylim=beta_range,bty='n')
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
  
}

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

# Split trials in each phase -----------------------------------------
# Aggregating behavior
Nalpha <- length(unique(Data_alpha$sub))
Nbeta <- length(unique(Data_beta$sub))
conf_group <- with(Data,aggregate(cj,by=list(phase_block,manip,sub,condition,cor),mean))
names(conf_group) <- c('phase_block','manip','sub','condition','cor','cj')

trials_phase <- data.frame(phase_block=rep(0:(Nphase_block-1),each=4),
                           cor=c(0,1),condition=rep(c("minus","plus"),each=2),sub=rep(subs,each=Nphase_block*4))

conf_group <- merge(conf_group,trials_phase,all=T)
conf_group[conf_group$sub %in% Data_alpha$sub,'manip'] <- 'alpha'
conf_group[conf_group$sub %in% Data_beta$sub,'manip'] <- 'beta'
table(complete.cases(conf_group$cj))
count_complete <- function(dat) {
  return(sum(complete.cases(dat)))
}

conf_alpha_count <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor,count_complete)
conf_alpha <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_alpha_sd <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_alpha_sd[,2:5] <- conf_alpha_sd[,2:5]/sqrt(conf_alpha_count[,2:5])
conf_beta <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_count <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,count_complete)
conf_beta_sd <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_beta_sd[,2:5] <- conf_beta_sd[,2:5]/sqrt(conf_beta_count[,2:5])

# Aggregating model fits
for (m in models) {
  print(m)
  cj_pred <- with(subset(anal_sim,model==m),aggregate(cj,by=list(trial,sub),mean))
  names(cj_pred) <- c("trial","sub","cj")
  cj_pred <- merge(cj_pred,Data[,c("trial","sub","phase_block","manip","condition","cor")])
  alpha <- with(subset(anal_sim,model==m),aggregate(alpha,by=list(trial,sub),mean))
  names(alpha) <- c("trial","sub","alpha")
  alpha <- merge(alpha,Data[,c("trial","sub","phase_block","manip","condition","cor")])
  beta <- with(subset(anal_sim,model==m),aggregate(beta,by=list(trial,sub),mean))
  names(beta) <- c("trial","sub","beta")
  beta <- merge(beta,Data[,c("trial","sub","phase_block","manip","condition","cor")])
  
  conf_group_pred <- with(cj_pred,aggregate(cj,by=list(phase_block,manip,sub,condition,cor),mean))
  names(conf_group_pred) <- c('phase_block','manip','sub','condition','cor','cj')
  
  conf_alpha_pred_count <- cast(subset(conf_group_pred,manip=='alpha'),phase_block~condition+cor, fun.aggregate = length)
  conf_alpha_pred <- cast(subset(conf_group_pred,manip=='alpha'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
  conf_alpha_pred_sd <- cast(subset(conf_group_pred,manip=='alpha'),phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
  conf_alpha_pred_sd[,2:5] <- conf_alpha_pred_sd[,2:5]/sqrt(conf_alpha_pred_count[,2:5])
  
  conf_beta_pred_count <- cast(subset(conf_group_pred,manip=='beta'),phase_block~condition+cor, fun.aggregate = length)
  conf_beta_pred <- cast(subset(conf_group_pred,manip=='beta'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
  conf_beta_pred_sd <- cast(subset(conf_group_pred,manip=='beta'),phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
  conf_beta_pred_sd[,2:5] <- conf_beta_pred_sd[,2:5]/sqrt(conf_beta_pred_count[,2:5])
  
  alpha_group_pred <- with(alpha,aggregate(alpha,by=list(phase_block,manip,sub,condition),mean))
  names(alpha_group_pred) <- c('phase_block','manip','sub','condition','cj')
  
  alpha_alpha_pred_count <- cast(subset(alpha_group_pred,manip=='alpha'),phase_block~condition, fun.aggregate = length)
  alpha_alpha_pred <- cast(subset(alpha_group_pred,manip=='alpha'),phase_block~condition,fun.aggregate = mean,na.rm=T)
  alpha_alpha_pred_sd <- cast(subset(alpha_group_pred,manip=='alpha'),phase_block~condition,fun.aggregate = sd,na.rm=T)
  alpha_alpha_pred_sd[,2:3] <- alpha_alpha_pred_sd[,2:3]/sqrt(alpha_alpha_pred_count[,2:3])
  
  alpha_beta_pred_count <- cast(subset(alpha_group_pred,manip=='beta'),phase_block~condition, fun.aggregate = length)
  alpha_beta_pred <- cast(subset(alpha_group_pred,manip=='beta'),phase_block~condition,fun.aggregate = mean,na.rm=T)
  alpha_beta_pred_sd <- cast(subset(alpha_group_pred,manip=='beta'),phase_block~condition,fun.aggregate = sd,na.rm=T)
  alpha_beta_pred_sd[,2:3] <- alpha_beta_pred_sd[,2:3]/sqrt(alpha_beta_pred_count[,2:3])
  
  
  beta_group_pred <- with(beta,aggregate(beta,by=list(phase_block,manip,sub,condition),mean))
  names(beta_group_pred) <- c('phase_block','manip','sub','condition','cj')
  
  beta_alpha_pred_count <- cast(subset(beta_group_pred,manip=='alpha'),phase_block~condition, fun.aggregate = length)
  beta_alpha_pred <- cast(subset(beta_group_pred,manip=='alpha'),phase_block~condition,fun.aggregate = mean,na.rm=T)
  beta_alpha_pred_sd <- cast(subset(beta_group_pred,manip=='alpha'),phase_block~condition,fun.aggregate = sd,na.rm=T)
  beta_alpha_pred_sd[,2:3] <- beta_alpha_pred_sd[,2:3]/sqrt(beta_alpha_pred_count[,2:3])
  
  beta_beta_pred_count <- cast(subset(beta_group_pred,manip=='beta'),phase_block~condition, fun.aggregate = length)
  beta_beta_pred <- cast(subset(beta_group_pred,manip=='beta'),phase_block~condition,fun.aggregate = mean,na.rm=T)
  beta_beta_pred_sd <- cast(subset(beta_group_pred,manip=='beta'),phase_block~condition,fun.aggregate = sd,na.rm=T)
  beta_beta_pred_sd[,2:3] <- beta_beta_pred_sd[,2:3]/sqrt(beta_beta_pred_count[,2:3])
  
  xlen <- nrow(conf_alpha)
  go_to("plots")
  tiff(paste0('trace_aggreg_alpha_',m,'_learn.tiff'), width = 36, height = 36, units = 'cm', res = 300)
  layout(matrix(c(1,2,1,3),ncol=2))
  par(mar=c(5,5,4,2)+.1)
  
  plot(conf_alpha$minus_0,ylim=c(.7,.9),col = BLUE, type = 'b',main=paste("Alpha-Manipulated Feedback",", empirical fits,",m,"learn"),
       lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",cex.main=cex.title,
       xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
  axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
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
  legend("top",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
         pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)
  
  par(mar=c(5,5,0,2)+.1)
  
  plot(alpha_alpha_pred$minus,ylim=c(5,15),col = BLUE, type = 'b',main="",
       lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Alpha",cex.main=cex.title,
       xlab = "",cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
  axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
  lines(alpha_alpha_pred$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
  error.bar(1:xlen,alpha_alpha_pred$minus,alpha_alpha_pred_sd$minus,
            lwd=2, col = BLUE)
  error.bar(1:xlen,alpha_alpha_pred$plus,alpha_alpha_pred_sd$plus,
            lwd=2, col = VERMILLION)
  
  plot(beta_alpha_pred$minus,ylim=c(5,15),col = BLUE, type = 'b',main="",
       lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Beta",cex.main=cex.title,
       xlab = "",cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
  axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
  lines(beta_alpha_pred$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
  error.bar(1:xlen,beta_alpha_pred$minus,beta_alpha_pred_sd$minus,
            lwd=2, col = BLUE)
  error.bar(1:xlen,beta_alpha_pred$plus,beta_alpha_pred_sd$plus,
            lwd=2, col = VERMILLION)
  
  dev.off()
  
  tiff(paste0('trace_aggreg_beta_',m,'_learn.tiff'), width = 36, height = 36, units = 'cm', res = 300)
  layout(matrix(c(1,2,1,3),ncol=2))
  par(mar=c(5,5,4,2)+.1)
  
  plot(conf_beta$minus_0,ylim=c(.7,.9),col = BLUE, type = 'b',main=paste("Beta-Manipulated Feedback",", empirical fits,",m,"learn"),
       lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",cex.main=cex.title,
       xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
  axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
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
  
  par(mar=c(5,5,0,2)+.1)
  
  plot(alpha_beta_pred$minus,ylim=c(5,15),col = BLUE, type = 'b',main="",
       lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Alpha",cex.main=cex.title,
       xlab = "",cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
  axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
  lines(alpha_beta_pred$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
  error.bar(1:xlen,alpha_beta_pred$minus,alpha_beta_pred_sd$minus,
            lwd=2, col = BLUE)
  error.bar(1:xlen,alpha_beta_pred$plus,alpha_beta_pred_sd$plus,
            lwd=2, col = VERMILLION)
  
  plot(beta_beta_pred$minus,ylim=c(5,15),col = BLUE, type = 'b',main="",
       lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Beta",cex.main=cex.title,
       xlab = "",cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
  axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
  lines(beta_beta_pred$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
  error.bar(1:xlen,beta_beta_pred$minus,beta_beta_pred_sd$minus,
            lwd=2, col = BLUE)
  error.bar(1:xlen,beta_beta_pred$plus,beta_beta_pred_sd$plus,
            lwd=2, col = VERMILLION)
  
  dev.off()
}


par(mfrow=c(1,1))


# Accuracy and RT over time -----------------------------------------------

# Aggregating behavior
count_complete <- function(dat) {
  return(sum(complete.cases(dat)))
}

if (is.factor(Data$cor)) {
  Data$cor <- as.numeric(Data$cor)-1
}

cex.legend <- 3

Nalpha <- length(unique(Data_alpha$sub))
Nbeta <- length(unique(Data_beta$sub))

trials_phase <- data.frame(phase_block=rep(0:(Nphase_block-1),each=2),
                           condition=c("minus","plus"),sub=rep(subs,each=Nphase_block*2))

rt_group <- with(Data,aggregate(rt,by=list(phase_block,manip,sub,condition),mean))
names(rt_group) <- c('phase_block','manip','sub','condition','cj')
rt_group <- merge(rt_group,trials_phase,all=T)
rt_group[rt_group$sub %in% Data_alpha$sub,'manip'] <- 'alpha'
rt_group[rt_group$sub %in% Data_beta$sub,'manip'] <- 'beta'

rt_alpha_count <- cast(subset(rt_group,manip=='alpha'),phase_block~condition,count_complete)
rt_alpha <- cast(subset(rt_group,manip=='alpha'),phase_block~condition,fun.aggregate = mean,na.rm=T)
rt_alpha_sd <- cast(subset(rt_group,manip=='alpha'),phase_block~condition,fun.aggregate = sd,na.rm=T)
rt_alpha_sd[,2:3] <- rt_alpha_sd[,2:3]/sqrt(rt_alpha_count[,2:3])
rt_beta <- cast(subset(rt_group,manip=='beta'),phase_block~condition,fun.aggregate = mean,na.rm=T)
rt_beta_count <- cast(subset(rt_group,manip=='beta'),phase_block~condition,count_complete)
rt_beta_sd <- cast(subset(rt_group,manip=='beta'),phase_block~condition,fun.aggregate = sd,na.rm=T)
rt_beta_sd[,2:3] <- rt_beta_sd[,2:3]/sqrt(rt_beta_count[,2:3])

cor_group <- with(Data,aggregate(cor,by=list(phase_block,manip,sub,condition),mean))
names(cor_group) <- c('phase_block','manip','sub','condition','cj')
cor_group <- merge(cor_group,trials_phase,all=T)
cor_group[cor_group$sub %in% Data_alpha$sub,'manip'] <- 'alpha'
cor_group[cor_group$sub %in% Data_beta$sub,'manip'] <- 'beta'

cor_alpha_count <- cast(subset(cor_group,manip=='alpha'),phase_block~condition,count_complete)
cor_alpha <- cast(subset(cor_group,manip=='alpha'),phase_block~condition,fun.aggregate = mean,na.rm=T)
cor_alpha_sd <- cast(subset(cor_group,manip=='alpha'),phase_block~condition,fun.aggregate = sd,na.rm=T)
cor_alpha_sd[,2:3] <- cor_alpha_sd[,2:3]/sqrt(cor_alpha_count[,2:3])
cor_beta <- cast(subset(cor_group,manip=='beta'),phase_block~condition,fun.aggregate = mean,na.rm=T)
cor_beta_count <- cast(subset(cor_group,manip=='beta'),phase_block~condition,count_complete)
cor_beta_sd <- cast(subset(cor_group,manip=='beta'),phase_block~condition,fun.aggregate = sd,na.rm=T)
cor_beta_sd[,2:3] <- cor_beta_sd[,2:3]/sqrt(cor_beta_count[,2:3])

xlen <- nrow(rt_alpha)
go_to("plots")
go_to("paper")
tiff('trace_rt_cor_.tiff', width = 36, height = 36, units = 'cm', res = 600)
layout(matrix(c(1,2,3,4),ncol=2,byrow = F))

plot(rt_alpha$minus,ylim=c(.8,1.2),col = BLUE, type = 'b',main="Alpha-Manipulated Feedback",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
title(ylab = "Reaction time (s)", xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"), 
      line = 2.5,cex.lab=cex.lab*.66/.83)
lines(rt_alpha$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
error.bar(1:xlen,rt_alpha$minus,rt_alpha_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:xlen,rt_alpha$plus,rt_alpha_sd$plus,
          lwd=2, col = VERMILLION)
# polygon(c(1:xlen,xlen:1),c(rt_alpha_pred$minus + rt_alpha_pred_sd$minus,
#                            (rt_alpha_pred$minus - rt_alpha_pred_sd$minus)[xlen:1]),
#         border=F,col=rgb(0,114,178,51,maxColorValue = 255))
# polygon(c(1:xlen,xlen:1),c(rt_alpha_pred$plus + rt_alpha_pred_sd$plus,
#                            (rt_alpha_pred$plus - rt_alpha_pred_sd$plus)[xlen:1]),
#         border=F,col=rgb(213,94,0,51,maxColorValue = 255))
legend("top",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)

plot(cor_alpha$minus,ylim=c(.5,1),col = BLUE, type = 'b',main="",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
title(ylab = "Accuracy", xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"), 
      line = 2.5,cex.lab=cex.lab*.66/.83)
lines(cor_alpha$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
error.bar(1:xlen,cor_alpha$minus,cor_alpha_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:xlen,cor_alpha$plus,cor_alpha_sd$plus,
          lwd=2, col = VERMILLION)
# polygon(c(1:xlen,xlen:1),c(cor_alpha_pred$minus + cor_alpha_pred_sd$minus,
#                            (cor_alpha_pred$minus - cor_alpha_pred_sd$minus)[xlen:1]),
#         border=F,col=rgb(0,114,178,51,maxColorValue = 255))
# polygon(c(1:xlen,xlen:1),c(cor_alpha_pred$plus + cor_alpha_pred_sd$plus,
#                            (cor_alpha_pred$plus - cor_alpha_pred_sd$plus)[xlen:1]),
#         border=F,col=rgb(213,94,0,51,maxColorValue = 255))

plot(rt_beta$minus,ylim=c(.8,1.2),col = BLUE, type = 'b',main="Beta-Manipulated Feedback",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
title(ylab = "Reaction time (s)", xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"), 
      line = 2.5,cex.lab=cex.lab*.66/.83)
lines(rt_beta$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
error.bar(1:xlen,rt_beta$minus,rt_beta_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:xlen,rt_beta$plus,rt_beta_sd$plus,
          lwd=2, col = VERMILLION)
# polygon(c(1:xlen,xlen:1),c(rt_beta_pred$minus + rt_beta_pred_sd$minus,
#                            (rt_beta_pred$minus - rt_beta_pred_sd$minus)[xlen:1]),
#         border=F,col=rgb(0,114,178,51,maxColorValue = 255))
# polygon(c(1:xlen,xlen:1),c(rt_beta_pred$plus + rt_beta_pred_sd$plus,
#                            (rt_beta_pred$plus - rt_beta_pred_sd$plus)[xlen:1]),
#         border=F,col=rgb(213,94,0,51,maxColorValue = 255))
legend("top",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)

plot(cor_alpha$minus,ylim=c(.5,1),col = BLUE, type = 'b',main="",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
title(ylab = "Accuracy", xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"), 
      line = 2.5,cex.lab=cex.lab*.66/.83)
lines(cor_alpha$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
error.bar(1:xlen,cor_alpha$minus,cor_alpha_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:xlen,cor_alpha$plus,cor_alpha_sd$plus,
          lwd=2, col = VERMILLION)
# polygon(c(1:xlen,xlen:1),c(cor_alpha_pred$minus + cor_alpha_pred_sd$minus,
#                            (cor_alpha_pred$minus - cor_alpha_pred_sd$minus)[xlen:1]),
#         border=F,col=rgb(0,114,178,51,maxColorValue = 255))
# polygon(c(1:xlen,xlen:1),c(cor_alpha_pred$plus + cor_alpha_pred_sd$plus,
#                            (cor_alpha_pred$plus - cor_alpha_pred_sd$plus)[xlen:1]),
#         border=F,col=rgb(213,94,0,51,maxColorValue = 255))
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



# Dynamic confidence plot according to group --------------------------------------------------------
m <- 'both'

conf_beta <- cast(subset(conf_group,manip=='beta' & sub%in%sub_nolearn_best),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_count <- cast(subset(conf_group,manip=='beta' & sub%in%sub_nolearn_best),phase_block~condition+cor,count_complete)
conf_beta_sd <- cast(subset(conf_group,manip=='beta' & sub%in%sub_nolearn_best),phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_beta_sd[,2:5] <- conf_beta_sd[,2:5]/sqrt(conf_beta_count[,2:5])

cj_pred <- with(subset(anal_sim,model==m & sub%in%sub_nolearn_best),aggregate(cj,by=list(trial,sub),mean))
names(cj_pred) <- c("trial","sub","cj")
cj_pred <- merge(cj_pred,Data[,c("trial","sub","phase_block","manip","condition","cor")])

conf_group_pred <- with(cj_pred,aggregate(cj,by=list(phase_block,manip,sub,condition,cor),mean))
names(conf_group_pred) <- c('phase_block','manip','sub','condition','cor','cj')

conf_beta_pred_count <- cast(subset(conf_group_pred,manip=='beta'),phase_block~condition+cor, fun.aggregate = length)
conf_beta_pred <- cast(subset(conf_group_pred,manip=='beta'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_pred_sd <- cast(subset(conf_group_pred,manip=='beta'),phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_beta_pred_sd[,2:5] <- conf_beta_pred_sd[,2:5]/sqrt(conf_beta_pred_count[,2:5])

cex_size <- function(size,cex.layout) {
  return(size/(par()$ps*cex.layout))
}
### Adjust sizes and positions
cex.layout <- 1
# Text
cex.legend <- cex_size(14,cex.layout)
cex.axis <- cex_size(14,cex.layout) 
cex.lab <- cex_size(18,cex.layout)
line.width <- 1
cex.pt <- .8


jpeg("conf_beta_bothlearn_nolearnsub.jpeg",width=13.5,height=10,units="cm",res=300)
par(mar=c(3,3,0,0)+.1)
plot(conf_beta$minus_0,ylim=c(.7,.915),col = BLUE, type = 'l',main="",
     lty = 2,  lwd = line.width, bty = 'n', xaxt = 'n', xlab="",ylab="",cex.axis=cex.axis)
title(ylab = "Confidence",xlab = "Within phase groups of 7 trials",cex.lab=cex.lab,line=2)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis)
lines(conf_beta$plus_0, type = 'l',  col = VERMILLION, lwd = line.width, lty = 2)
lines(conf_beta$plus_1, type = 'l',  col = VERMILLION, lwd = line.width, lty = 1)
lines(conf_beta$minus_1, type = 'l',  col = BLUE, lwd = line.width, lty = 1)
points(conf_beta$plus_0, pch = 16, cex = cex.pt, col = VERMILLION)
points(conf_beta$plus_1, pch = 16, cex = cex.pt, col = VERMILLION)
points(conf_beta$minus_1, pch = 16, cex = cex.pt, col = BLUE)
points(conf_beta$minus_0, pch = 16, cex = cex.pt, col = BLUE)
error.bar(1:length(conf_beta$minus_0),conf_beta$minus_0,conf_beta_sd$minus_0,
          lwd = line.width, col = BLUE)
error.bar(1:length(conf_beta$plus_0),conf_beta$plus_0,conf_beta_sd$plus_0,
          lwd = line.width, col = VERMILLION)
error.bar(1:length(conf_beta$minus_1),conf_beta$minus_1,conf_beta_sd$minus_1,
          lwd = line.width, col = BLUE)
error.bar(1:length(conf_beta$plus_1),conf_beta$plus_1,conf_beta_sd$plus_1,
          lwd = line.width, col = VERMILLION)
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
legend(0.5,.93,legend = c('High feedback','Low feedback'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16), bty = 'n',cex = cex.legend)
legend(12,.93,legend=c("Correct trials","Error trials"),
       title = NULL,lty=c(1,2),bty = "n",inset=0,
       cex = cex.legend, lwd=lwd.dat,seg.len=1.5)

dev.off()


# Dynamic confidence plot according to group 2 --------------------------------------------------------
m <- 'both'

conf_alpha <- cast(subset(conf_group,manip=='alpha' & sub%in%sub_learn_best),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_alpha_count <- cast(subset(conf_group,manip=='alpha' & sub%in%sub_learn_best),phase_block~condition+cor,count_complete)
conf_alpha_sd <- cast(subset(conf_group,manip=='alpha' & sub%in%sub_learn_best),phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_alpha_sd[,2:5] <- conf_alpha_sd[,2:5]/sqrt(conf_alpha_count[,2:5])

cj_pred <- with(subset(anal_sim,model==m & sub%in%sub_learn_best),aggregate(cj,by=list(trial,sub),mean))
names(cj_pred) <- c("trial","sub","cj")
cj_pred <- merge(cj_pred,Data[,c("trial","sub","phase_block","manip","condition","cor")])

conf_group_pred <- with(cj_pred,aggregate(cj,by=list(phase_block,manip,sub,condition,cor),mean))
names(conf_group_pred) <- c('phase_block','manip','sub','condition','cor','cj')

conf_alpha_pred_count <- cast(subset(conf_group_pred,manip=='alpha'),phase_block~condition+cor, fun.aggregate = length)
conf_alpha_pred <- cast(subset(conf_group_pred,manip=='alpha'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_alpha_pred_sd <- cast(subset(conf_group_pred,manip=='alpha'),phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_alpha_pred_sd[,2:5] <- conf_alpha_pred_sd[,2:5]/sqrt(conf_alpha_pred_count[,2:5])

cex_size <- function(size,cex.layout) {
  return(size/(par()$ps*cex.layout))
}
### Adjust sizes and positions
cex.layout <- 1
# Text
cex.legend <- cex_size(14,cex.layout)
cex.axis <- cex_size(14,cex.layout) 
cex.lab <- cex_size(18,cex.layout)
line.width <- 1
cex.pt <- .8


jpeg("conf_alpha_bothlearn_learnsub.jpeg",width=13.5,height=10,units="cm",res=300)
par(mar=c(3,3,0,0)+.1)
plot(conf_alpha$minus_0,ylim=c(.7,.915),col = BLUE, type = 'l',main="",
     lty = 2,  lwd = line.width, bty = 'n', xaxt = 'n', xlab="",ylab="",cex.axis=cex.axis)
title(ylab = "Confidence",xlab = "Within phase groups of 7 trials",cex.lab=cex.lab,line=2)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis)
lines(conf_alpha$plus_0, type = 'l',  col = VERMILLION, lwd = line.width, lty = 2)
lines(conf_alpha$plus_1, type = 'l',  col = VERMILLION, lwd = line.width, lty = 1)
lines(conf_alpha$minus_1, type = 'l',  col = BLUE, lwd = line.width, lty = 1)
points(conf_alpha$plus_0, pch = 16, cex = cex.pt, col = VERMILLION)
points(conf_alpha$plus_1, pch = 16, cex = cex.pt, col = VERMILLION)
points(conf_alpha$minus_1, pch = 16, cex = cex.pt, col = BLUE)
points(conf_alpha$minus_0, pch = 16, cex = cex.pt, col = BLUE)
error.bar(1:length(conf_alpha$minus_0),conf_alpha$minus_0,conf_alpha_sd$minus_0,
          lwd = line.width, col = BLUE)
error.bar(1:length(conf_alpha$plus_0),conf_alpha$plus_0,conf_alpha_sd$plus_0,
          lwd = line.width, col = VERMILLION)
error.bar(1:length(conf_alpha$minus_1),conf_alpha$minus_1,conf_alpha_sd$minus_1,
          lwd = line.width, col = BLUE)
error.bar(1:length(conf_alpha$plus_1),conf_alpha$plus_1,conf_alpha_sd$plus_1,
          lwd = line.width, col = VERMILLION)
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
legend(0.5,.93,legend = c('High feedback','Low feedback'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16), bty = 'n',cex = cex.legend)
legend(12,.93,legend=c("Correct trials","Error trials"),
       title = NULL,lty=c(1,2),bty = "n",inset=0,
       cex = cex.legend, lwd=lwd.dat,seg.len=1.5)

dev.off()


# Reanalyze confidence according to group -----------------------------------------------------------
if (stat_tests) {
  Data$group <- "static"
  Data[Data$sub %in% sub_learn_best,]$group <- "dynamic"
  Data_beta$group <- "static"
  Data_beta[Data_beta$sub %in% sub_learn_best,]$group <- "dynamic"
  Data_alpha$group <- "static"
  Data_alpha[Data_alpha$sub %in% sub_learn_best,]$group <- "dynamic"
  
  control <- lmerControl(optimizer = "bobyqa")
  glmercontrol <- glmerControl(optimizer = "bobyqa")
  
  # Set contrast coding
  options(contrasts=c("contr.sum","contr.poly"))
  
  
  # Alpha
  cj.cond.acc.interaction.alpha.0 <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (cor + condition + condition:cor|sub),
                                          data = Data_alpha, REML = F,control = control)
  cj.cond.acc.interaction.alpha <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + condition:withinphasetrial:group + (cor + condition + condition:cor|sub),
                                        data = Data_alpha, REML = F,control = control)
  anova(cj.cond.acc.interaction.alpha.0,cj.cond.acc.interaction.alpha)
  plot(resid(cj.cond.acc.interaction.alpha),Data_alpha$cj) #Linearity
  leveneTest(residuals(cj.cond.acc.interaction.alpha) ~ Data_alpha$cor*Data_alpha$condition*Data_alpha$difflevel) #Homogeneity of variance
  qqmath(cj.cond.acc.interaction.alpha) #Normality
  anova(cj.cond.acc.interaction.alpha) #Results  
  vif(cj.cond.acc.interaction.alpha)
  # post-hoc
  cj.static.group <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (cor + condition + condition:cor|sub),
                          data = subset(Data_alpha,group=="static"), REML = F,control = control)
  vif(cj.static.group)
  anova(cj.static.group)
  post_hoc <- emmeans(cj.static.group, ~ condition:withinphasetrial)
  pairs(post_hoc)
  post_hoc <- emmeans(cj.static.group, ~ condition:withinphasetrial | cor)
  pairs(post_hoc)
  cj.dynamic.group <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (cor + condition + condition:cor|sub),
                           data = subset(Data_alpha,group=="dynamic"), REML = F,control = control)
  vif(cj.dynamic.group)
  anova(cj.dynamic.group)
  post_hoc <- emmeans(cj.dynamic.group, ~ condition:withinphasetrial)
  pairs(post_hoc)
  post_hoc <- emmeans(cj.dynamic.group, ~ condition:withinphasetrial | cor)
  pairs(post_hoc)
  # Beta
  cj.cond.acc.interaction.beta.0 <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (cor + condition + condition:cor|sub),
                                         data = Data_beta, REML = F,control = control)
  cj.cond.acc.interaction.beta <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + condition:withinphasetrial:group + (cor + condition + condition:cor|sub),
                                       data = Data_beta, REML = F,control = control)
  anova(cj.cond.acc.interaction.beta.0,cj.cond.acc.interaction.beta)
  plot(resid(cj.cond.acc.interaction.beta),Data_beta$cj) #Linearity
  leveneTest(residuals(cj.cond.acc.interaction.beta) ~ Data_beta$cor*Data_beta$condition*Data_beta$difflevel) #Homogeneity of variance
  qqmath(cj.cond.acc.interaction.beta) #Normality
  anova(cj.cond.acc.interaction.beta) #Results  
  vif(cj.cond.acc.interaction.beta)
  # post-hoc
  cj.static.group <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (cor + condition + condition:cor|sub),
                          data = subset(Data_beta,group=="static"), REML = F,control = control)
  anova(cj.static.group)
  cj.dynamic.group <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (cor + condition + condition:cor|sub),
                           data = subset(Data_beta,group=="dynamic"), REML = F,control = control)
  anova(cj.dynamic.group)
}


# Questionnaire analysis --------------------------------------------------
library("ordinal")
questions$manip <- 'none'
questions[questions$sub %in% Data_beta$sub,'manip'] <- 'beta'
questions[questions$sub %in% Data_alpha$sub,'manip'] <- 'alpha'
table(questions$manip)
questions$group <- "static"
questions[questions$sub %in% sub_learn_best,]$group <- "dynamic"
# Recode Likert responses to numeric
questions[questions$questionResp=="['d']",]$questionResp <- 1
questions[questions$questionResp=="['f']",]$questionResp <- 2
questions[questions$questionResp=="['g']",]$questionResp <- 3
questions[questions$questionResp=="['h']",]$questionResp <- 4
questions[questions$questionResp=="['j']",]$questionResp <- 5
questions$group <- as.factor(questions$group)
questions$manip <- as.factor(questions$manip)

# 1. 'Did you notice anything unusual in this experiment ? If so, please specify.'
question1 <- subset(questions,questionID==1)

# Translate Q1 responses to English
# Install and load reticulate package
library(deeplr)

translated_text <- sapply(question1$questionResp, function(text) {
  translate2(text = text, source_lang = "NL", target_lang = "EN", auth_key = "0b2b0a60-b371-44d8-ba91-f7902122f0e9:fx")
})
question1$translation <- translated_text

# List of subjects that seem to have noticed the manipulation
aware_sub <- c(
  5, 7, 15, 16, 20, 21, 27, 28, 30, 34, 36, 37, 39, 41, 43, 45, 47, 49, 54, 56, 
  57, 59, 60, 77, 87, 89, 92, 95, 96, 97, 99, 106, 108, 109, 111, 112, 118, 128)
length(aware_sub)

Data$aware <- 0
Data[Data$sub %in% aware_sub,]$aware <- 1

Data_beta$aware <- 0
Data_beta[Data_beta$sub %in% aware_sub,]$aware <- 1

Data_alpha$aware <- 0
Data_alpha[Data_alpha$sub %in% aware_sub,]$aware <- 1

length_unique <- function(x) length(unique(x))
with(Data,aggregate(sub,by=list(aware,group),length_unique))
with(Data_beta,aggregate(sub,by=list(aware,group),length_unique))
with(Data_alpha,aggregate(sub,by=list(aware,group),length_unique))

aware.test <- glmer(aware ~ group + (1|sub), data=Data_beta, family="binomial", control = glmercontrol)
Anova(aware.test)

# 2. 'Indicate to what extent the feedback you received at the beginning of each task 
# helped you assess how well you performed.'
# Not at all     Slightly     Moderately     Very     Extremely
# D              F            G              H       J
question2 <- subset(questions,questionID==2)
question2$questionResp <- as.numeric(question2$questionResp)
with(question2,aggregate(questionResp,by=list(group,manip),table))
# Test difference in response between groups
question2$questionResp <- as.ordered(question2$questionResp)

ordinal.reg <- clm(data=question2,questionResp~group*manip)
summary(ordinal.reg)


# 3. 'Indicate to what extent you feel that the feedback you received closely 
# tracked your actual performance.'
question3 <- subset(questions,questionID==3)
question3$questionResp <- as.numeric(question3$questionResp)
with(question3,aggregate(questionResp,by=list(group,manip),table))
# Test difference in response between groups
question3$questionResp <- as.ordered(question3$questionResp)
ordinal.reg <- clm(data=question3,questionResp~group*manip)
summary(ordinal.reg)
ordinal.reg.beta <- clm(data=subset(question3,sub %in% Data_beta$sub),questionResp~group)
summary(ordinal.reg.beta)
ordinal.reg.alpha <- clm(data=subset(question3,sub %in% Data_alpha$sub),questionResp~group)
summary(ordinal.reg.alpha)

# Check participants understood the feedback
question4 <- subset(questions,questionID==5)

correct_response <- (grepl('3',question4$questionResp) | grepl('5',question4$questionResp)) & 
  !(grepl('1',question4$questionResp) | grepl('2',question4$questionResp) | grepl('6',question4$questionResp))
table(correct_response)

question4$fb_understood <- ifelse(correct_response,1,0)

fb_understood_sub <- subset(question4,fb_understood==1)$sub

Data$fb_understood <- 0
Data[Data$sub %in% fb_understood_sub,]$fb_understood <- 1

Data_beta$fb_understood <- 0
Data_beta[Data_beta$sub %in% fb_understood_sub,]$fb_understood <- 1

Data_alpha$fb_understood <- 0
Data_alpha[Data_alpha$sub %in% fb_understood_sub,]$fb_understood <- 1

with(Data_beta,aggregate(sub,by=list(fb_understood,group),length_unique))
with(Data_alpha,aggregate(sub,by=list(fb_understood,group),length_unique))

fb_understood.test <- glmer(fb_understood ~ group + (1|sub), data=Data, family="binomial", control = glmercontrol)
Anova(fb_understood.test)

fb_understood.test.beta <- glmer(fb_understood ~ group + (1|sub), data=Data_beta, family="binomial", control = glmercontrol)
Anova(fb_understood.test.beta)

fb_understood.test.alpha <- glmer(fb_understood ~ group + (1|sub), data=Data_alpha, family="binomial", control = glmercontrol)
Anova(fb_understood.test.alpha)

# Plot parameters per group -----------------------------------------------


par_beta <- subset(par,manip=='beta')
par_beta$group <- "static"
par_beta[par_beta$sub %in% sub_learn_best,]$group <- "dynamic"

with(par_beta,aggregate(eta_a,by=list(group=group,model=model),mean))
with(par_beta,aggregate(eta_b,by=list(group=group,model=model),mean))
with(par_beta,aggregate(a0,by=list(group=group,model=model),mean))
with(par_beta,aggregate(b0,by=list(group=group,model=model),mean))

par_alpha <- subset(par,manip=='alpha')
par_alpha$group <- "static"
par_alpha[par_alpha$sub %in% sub_learn_best,]$group <- "dynamic"

with(par_alpha,aggregate(eta_a,by=list(group=group,model=model),mean))
with(par_alpha,aggregate(eta_b,by=list(group=group,model=model),mean))
with(par_alpha,aggregate(a0,by=list(group=group,model=model),mean))
with(par_alpha,aggregate(b0,by=list(group=group,model=model),mean))

# Plot learning rates (eta_a and eta_b) separately for each group
par_beta$group <- factor(par_beta$group,levels=c("static","dynamic"))
par_alpha$group <- factor(par_alpha$group,levels=c("static","dynamic"))


ggplot(subset(par_beta,model=='both'), aes(x=eta_a, fill=group)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  # theme_ipsum() +
  labs(fill="")

ggplot(subset(par_beta,model=='both'), aes(x=eta_b, fill=group)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  # theme_ipsum() +
  labs(fill="")

ggplot(subset(par_alpha,model=='both'), aes(x=eta_a, fill=group)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  # theme_ipsum() +
  labs(fill="")

ggplot(subset(par_alpha,model=='both'), aes(x=eta_b, fill=group)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  # theme_ipsum() +
  labs(fill="")


# Learning rates
ggplot(subset(par_alpha,model=='both'), aes(x=eta_a, y=eta_b, color=group)) + 
  geom_point(size=2, alpha = .5) +
  ggtitle("Alpha experiment")

ggplot(subset(par_beta,model=='both'), aes(x=eta_a, y=eta_b, color=group)) +
  geom_point(size=2, alpha = .5) +
  ggtitle("Beta experiment")

# Initial values
ggplot(subset(par_alpha,model=='both'), aes(x=a0, y=b0, color=group)) + 
  geom_point(size=2, alpha = .5) +
  ggtitle("Alpha experiment")

ggplot(subset(par_beta,model=='both'), aes(x=a0, y=b0, color=group)) +
  geom_point(size=2, alpha = .5) +
  ggtitle("Beta experiment")


# Behavior results plot ---------------------------------------------------

# Behavior results figure layout ----------------------------------------------------
diff_order <- c('hard','medium','easy')
Ndiff <- length(diff_order)
Nbeta <- length(unique(Data_beta$sub))
cexax <- 1
cex.lab <- 1
cexleg <- 1
cex.title <- 1
cex.axis <- 1
linelab <- 2.5
lwdmean <- 3
mar.rt.acc <- c(2,3.5,0,1)
mar.cj <- c(4,3.5,0,1)
jpeg("Behavior_results.jpg", width = 19, height = 19, units = "cm", pointsize = 12, res = 1000)

layout(matrix(c(1,2,3,4,5,6),ncol=2,byrow = F),widths = c(1,2),heights = c(1,1,1.5))

par(mar=mar.rt.acc)

###
# RT
###

plot_rt_low <- with(subset(Data_beta,condition=="minus"&cor==1),
                    aggregate(rt,by=list(sub=sub,difflevel=difflevel),mean))
plot_rt_low <- cast(plot_rt_low,sub~difflevel)
plot_rt_low <- plot_rt_low[,diff_order] #Reorder columns to have hard -> easy

plot_rt_high <- with(subset(Data_beta,condition=="plus"&cor==1),
                     aggregate(rt,by=list(sub=sub,difflevel=difflevel),mean))
plot_rt_high <- cast(plot_rt_high,sub~difflevel)
plot_rt_high <- plot_rt_high[,diff_order] #Reorder columns to have hard -> easy

yrange <- c(min(c(colMeans(plot_rt_low),colMeans(plot_rt_high))) - max(c(colSds(plot_rt_low),colSds(plot_rt_high)))/sqrt(Nbeta),
            max(c(colMeans(plot_rt_low),colMeans(plot_rt_high))) + max(c(colSds(plot_rt_low),colSds(plot_rt_high)))/sqrt(Nbeta))
plot(xlab="",ylab="",colMeans(plot_rt_low),frame=F,type='n',xlim=c(.8,Ndiff+.2),
     ylim=yrange,xaxt='n')
title(ylab="Reaction time (s)",xlab="", line = linelab, cex.lab = cex.lab)
# axis(1,1:Ndiff,c("Hard","Average","Easy"),cex.axis=cexax)
axis(1,1:Ndiff,c("","",""),cex.axis=cexax)
points(colMeans(plot_rt_low),type='b',lwd=lwdmean,col=BLUE)
error.bar(1:Ndiff,colMeans(plot_rt_low),
          colSds(plot_rt_low,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=BLUE)

points(colMeans(plot_rt_high),type='b',lwd=lwdmean,col=VERMILLION)
error.bar(1:Ndiff,colMeans(plot_rt_high),
          colSds(plot_rt_high,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=VERMILLION)
legend("top",border=F,legend=c("Low","High"),pch=16,horiz=T,
       col=c(BLUE,VERMILLION),bty="n",
       cex=cexleg,title = "Feedback condition")
###
# Accuracy
###
par(mar=mar.rt.acc)


plot_cor_low <- with(subset(Data_beta,condition=="minus"),
                     aggregate(cor,by=list(sub=sub,difflevel=difflevel),mean))
plot_cor_low <- cast(plot_cor_low,sub~difflevel)
plot_cor_low <- plot_cor_low[,diff_order] #Reorder columns to have hard -> easy

plot_cor_high <- with(subset(Data_beta,condition=="plus"),
                      aggregate(cor,by=list(sub=sub,difflevel=difflevel),mean))
plot_cor_high <- cast(plot_cor_high,sub~difflevel)
plot_cor_high <- plot_cor_high[,diff_order] #Reorder columns to have hard -> easy

yrange <- c(min(colMeans(rbind(plot_cor_low,plot_cor_high))) - min(colSds(rbind(plot_cor_low,plot_cor_high)))/sqrt(Nbeta),
            max(colMeans(rbind(plot_cor_low,plot_cor_high))) + max(colSds(rbind(plot_cor_low,plot_cor_high)))/sqrt(Nbeta))
plot(xlab="",ylab="",colMeans(plot_cor_low),frame=F,type='n',xlim=c(.8,Ndiff+.2),
     ylim=c(.5,1),xaxt='n')
title(ylab="Accuracy",xlab="", line = linelab, cex.lab = cex.lab)
axis(1,1:Ndiff,c("","",""),cex.axis=cexax)
points(colMeans(plot_cor_low),type='b',lwd=lwdmean,col=BLUE)
error.bar(1:Ndiff,colMeans(plot_cor_low),
          colSds(plot_cor_low,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=BLUE)

points(colMeans(plot_cor_high),type='b',lwd=lwdmean,col=VERMILLION)
error.bar(1:Ndiff,colMeans(plot_cor_high),
          colSds(plot_cor_high,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=VERMILLION)
legend("top",border=F,legend=c("Low","High"),pch=16,horiz=T,
       col=c(BLUE,VERMILLION),bty="n",
       cex=cexleg,title = "Feedback condition")

###
# Confidence
###
par(mar=mar.cj)

plot_cj_low <- with(subset(Data_beta,condition=="minus"),
                    aggregate(cj_integer,by=list(sub=sub,difflevel=difflevel,cor=cor),mean))
plot_cj_low_cor <- cast(subset(plot_cj_low,cor==1),sub~difflevel)
plot_cj_low_err <- cast(subset(plot_cj_low,cor==0),sub~difflevel)
plot_cj_low_cor <- plot_cj_low_cor[,diff_order] #Reorder columns to have hard -> easy
plot_cj_low_err <- plot_cj_low_err[,diff_order] #Reorder columns to have hard -> easy

plot_cj_high <- with(subset(Data_beta,condition=="plus"),
                     aggregate(cj_integer,by=list(sub=sub,difflevel=difflevel,cor=cor),mean))
plot_cj_high_cor <- cast(subset(plot_cj_high,cor==1),sub~difflevel)
plot_cj_high_err <- cast(subset(plot_cj_high,cor==0),sub~difflevel)
plot_cj_high_cor <- plot_cj_high_cor[,diff_order] #Reorder columns to have hard -> easy
plot_cj_high_err <- plot_cj_high_err[,diff_order] #Reorder columns to have hard -> easy

yrange <- c(min(colMeans(rbind(plot_cj_low_err,plot_cj_high_err),na.rm = T)) - min(colSds(rbind(plot_cj_low_err,plot_cj_high_err),na.rm = T))/sqrt(Nbeta),
            max(colMeans(rbind(plot_cj_low_cor,plot_cj_high_cor),na.rm = T)) + max(colSds(rbind(plot_cj_low_cor,plot_cj_high_cor),na.rm = T))/sqrt(Nbeta))
plot(xlab="",ylab="",colMeans(plot_cj_low_cor),frame=F,type='n',xlim=c(.8,Ndiff+.2),
     ylim=yrange,xaxt='n')
title(ylab="Confidence",xlab="Trial difficulty", line = linelab, cex.lab = cex.lab)
axis(1,1:Ndiff,c("Hard","Average","Easy"),cex.axis=cexax)
# Incorrect trials
points(colMeans(plot_cj_low_err,na.rm=T),type='b',lwd=lwdmean,col=BLUE,lty=2)
error.bar(1:Ndiff,colMeans(plot_cj_low_err),
          colSds(plot_cj_low_err,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=BLUE,lty=1)
points(colMeans(plot_cj_high_err,na.rm=T),type='b',lwd=lwdmean,col=VERMILLION,lty=2)
error.bar(1:Ndiff,colMeans(plot_cj_high_err,na.rm=T),
          colSds(plot_cj_high_err,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=VERMILLION,lty=1)
# Correct trials
points(colMeans(plot_cj_low_cor,na.rm=T),type='b',lwd=lwdmean,col=BLUE,lty=1)
error.bar(1:Ndiff,colMeans(plot_cj_low_cor),
          colSds(plot_cj_low_cor,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=BLUE,lty=1)
points(colMeans(plot_cj_high_cor,na.rm=T),type='b',lwd=lwdmean,col=VERMILLION,lty=1)
error.bar(1:Ndiff,colMeans(plot_cj_high_cor,na.rm=T),
          colSds(plot_cj_high_cor,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=VERMILLION,lty=1)
legend("top",border=F,legend=c("Low","High"),pch=16,horiz=T,
       col=c(BLUE,VERMILLION),bty="n",
       cex=cexleg,title = "Feedback condition")

# Plot aggregated traces  -------------------------------------------------
m <- 'both'
cj_pred <- with(subset(anal_sim,model==m),aggregate(cj,by=list(trial,sub),mean))
names(cj_pred) <- c("trial","sub","cj")
cj_pred <- merge(cj_pred,Data[,c("trial","sub","phase_block","manip","condition","cor")])

conf_group_pred <- with(cj_pred,aggregate(cj,by=list(phase_block,manip,sub,condition,cor),mean))
names(conf_group_pred) <- c('phase_block','manip','sub','condition','cor','cj')

conf_beta_pred_count <- cast(subset(conf_group_pred,manip=='beta'),phase_block~condition+cor, fun.aggregate = length)
conf_beta_pred <- cast(subset(conf_group_pred,manip=='beta'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_pred_sd <- cast(subset(conf_group_pred,manip=='beta'),phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_beta_pred_sd[,2:5] <- conf_beta_pred_sd[,2:5]/sqrt(conf_beta_pred_count[,2:5])


if (is.factor(Data$cor)) {
  Data$cor <- as.numeric(Data$cor)-1
}


trials_phase <- data.frame(phase_block=rep(0:(Nphase_block-1),each=2),
                           condition=c("minus","plus"),sub=rep(subs,each=Nphase_block*2))

rt_group <- with(Data,aggregate(rt,by=list(phase_block,manip,sub,condition),mean))
names(rt_group) <- c('phase_block','manip','sub','condition','cj')
rt_group <- merge(rt_group,trials_phase,all=T)
rt_group[rt_group$sub %in% Data_beta$sub,'manip'] <- 'beta'

rt_beta <- cast(subset(rt_group,manip=='beta'),phase_block~condition,fun.aggregate = mean,na.rm=T)
rt_beta_count <- cast(subset(rt_group,manip=='beta'),phase_block~condition,fun.aggregate=length)
rt_beta_sd <- cast(subset(rt_group,manip=='beta'),phase_block~condition,fun.aggregate = sd,na.rm=T)
rt_beta_sd[,2:3] <- rt_beta_sd[,2:3]/sqrt(rt_beta_count[,2:3])

cor_group <- with(Data,aggregate(cor,by=list(phase_block,manip,sub,condition),mean))
names(cor_group) <- c('phase_block','manip','sub','condition','cj')
cor_group <- merge(cor_group,trials_phase,all=T)
cor_group[cor_group$sub %in% Data_beta$sub,'manip'] <- 'beta'


cor_beta <- cast(subset(cor_group,manip=='beta'),phase_block~condition,fun.aggregate = mean,na.rm=T)
cor_beta_count <- cast(subset(cor_group,manip=='beta'),phase_block~condition,fun.aggregate=length)
cor_beta_sd <- cast(subset(cor_group,manip=='beta'),phase_block~condition,fun.aggregate = sd,na.rm=T)
cor_beta_sd[,2:3] <- cor_beta_sd[,2:3]/sqrt(cor_beta_count[,2:3])

conf_group <- with(Data,aggregate(cj,by=list(phase_block,manip,sub,condition,cor),mean))
names(conf_group) <- c('phase_block','manip','sub','condition','cor','cj')
conf_group <- merge(conf_group,trials_phase,all=T)
conf_group[conf_group$sub %in% Data_alpha$sub,'manip'] <- 'alpha'
conf_group[conf_group$sub %in% Data_beta$sub,'manip'] <- 'beta'
conf_beta <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_count <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,length)
conf_beta_sd <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_beta_sd[,2:5] <- conf_beta_sd[,2:5]/sqrt(conf_beta_count[,2:5])


par(mar=mar.rt.acc)

xlen <- dim(conf_beta)[1]
plot(rt_beta$minus,ylim=c(.8,1.2),col = BLUE, type = 'b',main="",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab,cex.axis=cex.axis)
axis(1, at = 1:Nphase_block, labels = F,cex.axis=cex.axis)
title(ylab = "Reaction time (s)", xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"), 
      line = linelab,cex.lab=cex.lab)
lines(rt_beta$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
error.bar(1:xlen,rt_beta$minus,rt_beta_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:xlen,rt_beta$plus,rt_beta_sd$plus,
          lwd=2, col = VERMILLION)

par(mar=mar.rt.acc)

plot(cor_beta$minus,ylim=c(.5,1),col = BLUE, type = 'b',main="",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab,cex.axis=cex.axis)
axis(1, at = 1:Nphase_block, labels = F,cex.axis=cex.axis)
title(ylab = "Accuracy", xlab = "", 
      line = linelab,cex.lab=cex.lab)
lines(cor_beta$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
error.bar(1:xlen,cor_beta$minus,cor_beta_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:xlen,cor_beta$plus,cor_beta_sd$plus,
          lwd=2, col = VERMILLION)

par(mar=mar.cj)


plot(conf_beta$minus_0,ylim=c(.7,.9),col = BLUE, type = 'b',main="",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab,cex.axis=cex.axis)
title(ylab="Confidence",xlab= paste("Within phase groups of",Nphase_trial/Nphase_block,"trials"), line = linelab, cex.lab = cex.lab)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis)
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
legend("top",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n',cex = cexleg)

dev.off()
par(mar=c(5,4,4,2)+.1)

