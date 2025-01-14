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
sd(age_alpha$x)
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

# Let's check how many participants were excluded from chance performance
table(unique(exclusion) %in% subset(Data,manip=='beta')$sub)

#' Filter out participants who reported only one confidence level 
#' more than 90% of the time
conf_count <- with(Data,aggregate(cj,by=list(sub=sub),max_count))
conf_count$x <- conf_count$x/Ntrials
exclusion <- c(exclusion, unique(conf_count[conf_count$x>.9,"sub"]))

# Let's check how many participants were excluded from lack of confidence variability
table(unique(exclusion) %in% subset(Data,manip=='beta')$sub)

Data <- subset(Data,!(sub %in% exclusion))
questions <- subset(questions,!(sub %in% exclusion))

# Convert RTs to seconds
Data$rt <- Data$rt/1000
Data$RTconf <- Data$RTconf/1000

Data$response[Data$response==0] <- -1

Data$condition <- as.factor(Data$condition)
Data$difflevel <- as.factor(Data$difflevel)
subs <- unique(Data$sub); Nsub <- length(subs)
conditions <- sort(unique(Data$condition))
difflevels <- sort(unique(Data$difflevel))

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
Data <- subset(Data,RTconf<5 & rt<5 & rt>.2)
write.csv(Data,file = "alternating_fb_aggregated_data.csv",row.names = F)

Nphase_trial <- length(unique(Data$withinphasetrial))
Nphase_block <- 18
Data$phase_block <-  Data$withinphasetrial %/% (Nphase_trial/Nphase_block)
Data$phase_block <- as.factor(Data$phase_block)

Data$cj_integer <- Data$cj*6

Data_alpha <- subset(Data,manip=='alpha')
Data_beta <- subset(Data,manip=='beta')
subs_alpha <- unique(Data_alpha$sub); Nalpha <- length(subs_alpha)
subs_beta <- unique(Data_beta$sub); Nbeta <- length(subs_beta)

# Behavior analysis - Static ----------------------------------------------
if (stat_tests) {
  ##' We first look at the overall effect of feedback on behavior to see if we 
  ##' replicate the findings of previous studies.
  
  control <- lmerControl(optimizer = "bobyqa")
  glmercontrol <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5))
  
  # Set contrast coding
  options(contrasts=c("contr.sum","contr.poly"))
  
  ## Reaction Times 
  # Alpha experiment
  rt.int.alpha.static <- lmer(rt~condition*difflevel + (1|sub),data = subset(Data_alpha,cor==1),REML = F,control = control)
  rt.cond.alpha.static <- lmer(rt~condition*difflevel + (condition|sub),data = subset(Data_alpha,cor==1),REML = F,control = control)
  anova(rt.int.alpha.static,rt.cond.alpha.static)
  # Singular fit
  # rt.cond.diff.alpha.static <- lmer(rt~condition*difflevel + (condition+difflevel|sub),data = subset(Data_alpha,cor==1),REML = F,control = control)
  anova(rt.cond.alpha.static)
  
  # Beta experiment
  rt.int.beta.static <- lmer(rt~condition*difflevel + (1|sub),data = subset(Data_beta,cor==1),REML = F,control = control)
  rt.cond.beta.static <- lmer(rt~condition*difflevel + (condition|sub),data = subset(Data_beta,cor==1),REML = F,control = control)
  anova(rt.int.beta.static,rt.cond.beta.static)
  # Singular fit
  # rt.cond.diff.beta.static <- lmer(rt~condition*difflevel + (condition+difflevel|sub),data = subset(Data_beta,cor==1),REML = F,control = control)
  anova(rt.cond.beta.static)
  
  ## Accuracy
  # Alpha experiment
  acc.int.alpha.static <- glmer(cor~condition*difflevel + (1|sub),data = Data_alpha,family = binomial, control = glmercontrol)
  acc.cond.alpha.static <- glmer(cor~condition*difflevel + (condition|sub),data = Data_alpha,family = binomial, control = glmercontrol)
  anova(acc.int.alpha.static,acc.cond.alpha.static)
  acc.cond.diff.alpha.static <- glmer(cor~condition*difflevel + (condition+difflevel|sub),data = Data_alpha,family = binomial)
  Anova(acc.cond.alpha.static)
  
  # Beta experiment
  acc.int.beta.static <- glmer(cor~condition*difflevel + (1|sub),data = Data_beta,family = binomial, control = glmercontrol)
  acc.cond.beta.static <- glmer(cor~condition*difflevel + (condition|sub),data = Data_beta,family = binomial, control = glmercontrol)
  anova(acc.int.beta.static,acc.cond.beta.static)
  acc.cond.diff.beta.static <- glmer(cor~condition*difflevel + (condition+difflevel|sub),data = Data_beta,family = binomial)
  Anova(acc.cond.beta.static)
  
  ## Confidence
  # Alpha experiment
  cj.int.alpha.static <- lmer(cj ~ condition*cor*difflevel + (1|sub),data = Data_alpha,REML = F,control = control); 
  cj.cond.alpha.static <- lmer(cj ~ condition*cor*difflevel + (condition|sub),data = Data_alpha,REML = F,control = control); 
  anova(m.int,m.cond)
  cj.cond.acc.alpha.static <- lmer(cj ~ condition*cor*difflevel + (cor + condition|sub),data = Data_alpha, REML = F,control = control)
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
  
  # Beta experiment
  cj.int.beta.static <- lmer(cj ~ condition*cor*difflevel + (1|sub),data = Data_beta,REML = F,control = control); 
  cj.cond.beta.static <- lmer(cj ~ condition*cor*difflevel + (condition|sub),data = Data_beta,REML = F,control = control); 
  anova(cj.int.beta.static,cj.cond.beta.static)
  cj.cond.acc.beta.static <- lmer(cj ~ condition*cor*difflevel + (cor + condition|sub),data = Data_beta, REML = F,control = control)
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

# Behavior analysis - Dynamics -----------------------------------------------------------
if (stat_tests) {

  ## Reaction Times
  # Alpha experiment
  rt.int.alpha <- lmer(rt~condition*difflevel*withinphasetrial + (1|sub),data = subset(Data_alpha,cor==1),REML = F,control = control)
  rt.cond.alpha <- lmer(rt~condition*difflevel*withinphasetrial + (condition|sub),data = subset(Data_alpha,cor==1),REML = F,control = control)
  anova(rt.int.alpha,rt.cond.alpha)
  # Singular fit
  # rt.cond.diff.alpha <- lmer(rt~condition*difflevel*withinphasetrial + (condition+difflevel|sub),data = subset(Data_alpha,cor==1),REML = F,control = control)
  anova(rt.cond.alpha)
  
  # Beta experiment
  rt.int.beta <- lmer(rt~condition*difflevel*withinphasetrial + (1|sub),data = subset(Data_beta,cor==1),REML = F,control = control)
  rt.cond.beta <- lmer(rt~condition*difflevel*withinphasetrial + (condition|sub),data = subset(Data_beta,cor==1),REML = F,control = control)
  anova(rt.int.beta,rt.cond.beta)
  # Singular fit
  # rt.cond.diff.beta <- lmer(rt~condition*difflevel*withinphasetrial + (condition+difflevel|sub),data = subset(Data_beta,cor==1),REML = F,control = control)
  anova(rt.cond.beta)
  
  ## Accuracy
  # Alpha experiment
  acc.int.alpha <- glmer(cor~condition*difflevel*withinphasetrial + (1|sub),data = Data_alpha,family = binomial, control = glmercontrol)
  # Somehow even the simplest mixed model doesn't converge. Using regular Anova instead.
  acc.alpha <- glm(cor~condition*difflevel*withinphasetrial,data = Data_alpha,family = binomial)
  Anova(acc.alpha,type=3)
  
  # Beta experiment
  acc.int.beta <- glmer(cor~condition*difflevel*withinphasetrial + (1|sub),data = Data_beta,family = binomial, control = glmercontrol)
  # Same thing happened here. Using regular Anova again
  acc.beta <- glm(cor~condition*difflevel*withinphasetrial,data = Data_beta,family = binomial)
  Anova(acc.beta,type=3)
  
  
  ## Confidence
  # Alpha experiment
  cj.cond.interaction.alpha.cor <- lmer(cj ~ condition*difflevel*withinphasetrial + (condition|sub),
                                        data = subset(Data_alpha,cor==1), REML = F,control = control)
  anova(cj.cond.interaction.alpha.cor)

  
  # Beta experiment
  cj.int.beta <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (1|sub),data = Data_beta,REML = F,control = control); 
  cj.cond.beta <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (condition|sub),data = Data_beta,REML = F,control = control); 
  anova(cj.int.beta,cj.cond.beta)
  cj.cond.acc.beta <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + (cor + condition|sub),data = Data_beta, REML = F,control = control)
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


# Retrieve fits -----------------------------------------------------------
# Some fixed arguments
beta_input <- .1 # Beta scaled to .1 improves convergence
dt <- .001; sigma <- .1
binning <- F
Nsim <- 20

models <- c("no","alpha","beta","both")

totlen <- length(difflevels)*length(subs_beta)*length(models)

# DDM parameters + DLDC parameters + MSE
par <- data.frame(cost_ldc = NA, a0 = NA, b0 = NA, eta_a = NA, eta_b = NA,
                  bound = NA, drift = NA, ter = NA,cost_ddm = NA, 
                  sub = rep(subs_beta, each = length(difflevels)*length(models)),
                  difflevel = rep(difflevels, length.out = totlen),
                  manip = NA, model = rep(models, each = length(difflevels)) )

Data_beta$drift <- NA

Data_beta <- Data_beta[order(Data_beta$sub,Data_beta$trial),]
Data <- Data[order(Data$sub,Data$trial),]

go_to("fit")

# Create a data frame to store model simulations
if (!file.exists("anal_sim.Rdata")) {
  anal_sim <- Data_beta[rep(seq_len(nrow(Data_beta)), each=Nsim*length(models)), c('trial','withinphasetrial','sub','condition','cor')]
  anal_sim$alpha <- NA
  anal_sim$beta <- NA
  anal_sim$sim <- 1:Nsim
  anal_sim$model <- rep(models,each=Nsim)
}

# Retrieve model fits
for (s in 1:Nbeta) {
  
  print(paste("Retrieving participant",s))
  
  temp_dat <- subset(Data_beta,sub==subs_beta[s])
  
  # DDM fits
  ddm_file <- paste0('ddm/ddmfit_',subs_beta[s],'.Rdata')
  if(file.exists(ddm_file)){
    load(ddm_file)
  }else{
    next
  }
  ddm_params <- ddm.results$optim$bestmem[c(1,5:length(ddm.results$optim$bestmem))]
  
  # Store estimated parameters in big dataframe
  par[par$sub==subs_beta[s],"bound"] <- ddm.results$optim$bestmem[1]
  par[par$sub==subs_beta[s],"ter"] <- ddm.results$optim$bestmem[2]
  par[par$sub==subs_beta[s],"z"] <- ddm.results$optim$bestmem[3]
  par[par$sub==subs_beta[s],"vratio"] <- ddm.results$optim$bestmem[4]
  par[par$sub==subs_beta[s]&par$difflevel==difflevels[1],"drift"] <- ddm.results$optim$bestmem[5]
  par[par$sub==subs_beta[s]&par$difflevel==difflevels[2],"drift"] <- ddm.results$optim$bestmem[6]
  par[par$sub==subs_beta[s]&par$difflevel==difflevels[3],"drift"] <- ddm.results$optim$bestmem[7]
  par[par$sub==subs_beta[s],"cost_ddm"] <- ddm.results$optim$bestval
  
  # Add estimated drift to behavior data DF
  for (diff in difflevels) {
    Data_beta$drift[Data_beta$difflevel==diff&Data_beta$sub==subs_beta[s]] <- unique(par[par$sub==subs_beta[s]&par$difflevel==diff,"drift"])
  }
  
  temp_dat <- subset(Data_beta,sub==subs_beta[s])
  
  # Retrieve DLDC fits of all 4 variant models
  for (model in models) {
    ldc_file <- paste0('ldc_nn/behavior_fit/',model,'_learn/ldcfit_',subs_beta[s],'.Rdata')
    if (file.exists(ldc_file)) {
      load(ldc_file)
      par[par$sub==subs_beta[s]&par$model==model,"a0"] <- ldc.results$optim$bestmem[1]
      par[par$sub==subs_beta[s]&par$model==model,"b0"] <- ldc.results$optim$bestmem[2]
      par[par$sub==subs_beta[s]&par$model==model,"eta_a"] <- ldc.results$optim$bestmem[4]
      par[par$sub==subs_beta[s]&par$model==model,"eta_b"] <- ldc.results$optim$bestmem[5]
      par[par$sub==subs_beta[s]&par$model==model,"cost_ldc"] <- ldc.results$optim$bestval
      
      #' Generate model predictions multiple times to account for the stochastic
      #' nature of estimating single trial accumulated evidence 
      if (!file.exists("anal_sim.Rdata")) {
        for (i in 1:Nsim) {
          results <-
            ldc.nn.fit.w(params=c(mean(par[par$sub==subs_beta[s]&par$model==model,"a0"]),
                                  mean(par[par$sub==subs_beta[s]&par$model==model,"b0"]),1,
                                  mean(par[par$sub==subs_beta[s]&par$model==model,"eta_a"]),
                                  mean(par[par$sub==subs_beta[s]&par$model==model,"eta_b"])),
                         ddm_params = ddm_params,
                         obs=temp_dat,returnFit = F,eta_sep=T,
                         binning = binning,
                         dt = dt, sigma = sigma)
          anal_sim[anal_sim$sim==i&anal_sim$sub==subs_beta[s]&anal_sim$model==model ,'cj'] <- results$pred
          anal_sim[anal_sim$sim==i&anal_sim$sub==subs_beta[s]&anal_sim$model==model ,'alpha'] <- results$trace[,1]
          anal_sim[anal_sim$sim==i&anal_sim$sub==subs_beta[s]&anal_sim$model==model ,'beta'] <- results$trace[,2]  
        }
      }
    }
  }
  
  par[par$sub==subs_beta[s],"manip"] <- unique(temp_dat$manip)
}

# Save model simulations
if (!file.exists("anal_sim.Rdata")) {
  save(anal_sim,file="anal_sim.Rdata")
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
  Data_beta <- merge(Data_beta,anal_sim_mean,by=c("trial","sub"))
  
}
Data_beta <- Data_beta[order(Data_beta$sub,Data_beta$trial),]

anal_sim$cj_scaled <- anal_sim$cj*6

# Model comparison --------------------------------------------------------
par$Ndata_point <-  round(nrow(Data_beta)/Nbeta)
par$Npar <- 3
par[par$model=="no",'Npar'] <- 2
par[par$model=="both",'Npar'] <- 4

bic_custom <- function(Residuals,k,n){
  return(log(n)*k+n*log(Residuals/n))
}

### Mean BIC over participants per model
par$bic <- bic_custom(par$cost_ldc,par$Npar,par$Ndata_point)
mean_bic <- with(par,aggregate(bic,by=list(model=model),mean))
mean_bic$delta <- -99
mean_bic$delta <- mean_bic$x - min(mean_bic$x)


### Best model per participant
model_bic <- with(par,aggregate(bic,list(sub=sub,model=model),mean))
model_bic <- cast(model_bic, sub ~ model)
models_ord <- names(model_bic)[2:5]
model_bic$best <- models_ord[apply(model_bic[,2:5],1,which.min)]
model_bic$best_learning <- apply(model_bic[,3:5],1,min)
table(model_bic$best)
model_bic$diff_learn <- model_bic$best_learning - model_bic$no


# Plot differences in BIC
go_to("plots")
jpeg("model_comparison.jpg",width=19,height=6,units = 'cm',res=600, pointsize = 10)
par(mfrow=c(1,2))
par(mar=c(4,5.5,2,0))
mean_bic_nolearn <- subset(mean_bic,model=="no" )$delta
mean_bic_learn <- min(subset(mean_bic,model!="no" )$delta) 
range_bic_learn <- range(subset(mean_bic,model!="no" )$delta) 
bp <- barplot(matrix(c(mean_bic_learn,mean_bic_nolearn)), las = 2, horiz=T,xaxt='n',
              names.arg = c("Learning", "No learning"), border = NA, beside = T,
              xlab="",ylab = "")
title(xlab = expression(paste(Delta,"BIC")), ylab = "", line = 2.5)
axis(1, seq(0,14,2), seq(0,14,2))

# Plot N best models
par(mar=c(4,3.5,2,2))
if ("no" %in% model_bic$best) {
  model_bic[model_bic$best != 'no',"best"] <- "Learning"
  model_bic[model_bic$best == 'no',"best"] <- "No learning"
}
best_model <- table(model_bic$best)
bp <- barplot(matrix(best_model), las=2, xlab="",ylab = "",horiz=T, 
              xaxt='n',xlim = c(0,35), beside =T, names.arg = c("", ""), border = NA)
title(xlab = "Participant count", ylab = "", line = 2.5)
axis(1, seq(0,35,5), seq(0,35,5))
# Add model names under each bar
dev.off()
par(mar=c(5,4,4,2)+.1)
par(mfrow=c(1,1))

sub_nolearn_best <- subset(model_bic,best=='No learning')$sub
sub_learn_best <- subset(model_bic,best!='No learning')$sub

# Add the best model to the empirical data frame
if ("best" %in% names(Data_beta)) {
  Data_beta <- Data_beta[,!(names(Data_beta) == "best")]
}
Data_beta <- merge(Data_beta,model_bic[,c('sub','best','diff_learn')])
Data_beta$best <- as.factor(Data_beta$best)

# Create groups according to best model
Data_beta$group <- "dynamic"
Data_beta[Data_beta$best == "No learning",]$group <- "static"
# Analysis of model prediction --------------------------------------------
if (stat_tests) {
  
  cj.pred.cond.acc.interaction.beta.static <- lmer(cj_pred_both_learn ~ condition*cor*difflevel + (cor + condition + condition:cor|sub),
                                                   data = Data_beta, REML = F,control = control)
  plot(resid(cj.pred.cond.acc.interaction.beta.static),Data_beta$cj_pred_both_learn) #Linearity
  leveneTest(residuals(cj.pred.cond.acc.interaction.beta.static) ~ Data_beta$cor*Data_beta$condition*Data_beta$difflevel) #Homogeneity of variance
  qqmath(cj.pred.cond.acc.interaction.beta.static) #Normality
  anova(cj.pred.cond.acc.interaction.beta.static) #Results  
  post_hoc <- emmeans(cj.pred.cond.acc.interaction.beta.static, ~ difflevel)
  pairs(post_hoc)
  
  cj.pred.cond.acc.interaction.beta <- lmer(cj_pred_both_learn~ condition*cor*difflevel*withinphasetrial + (cor + condition + condition:cor|sub),
                                            data = Data_beta, REML = F,control = control)
  plot(resid(cj.pred.cond.acc.interaction.beta),Data_beta$cj_pred_both_learn) #Linearity
  leveneTest(residuals(cj.pred.cond.acc.interaction.beta) ~ Data_beta$cor*Data_beta$condition*Data_beta$difflevel) #Homogeneity of variance
  qqmath(cj.pred.cond.acc.interaction.beta) #Normality
  vif(cj.pred.cond.acc.interaction.beta) # Multicollinearity
  anova(cj.pred.cond.acc.interaction.beta) #Results  
}
# Compute confidence rolling mean per subject ----------------------------------------
n <- 25 # Rolling mean window size
n_err <- 25

trial_conf_sub <- with(Data,aggregate(cj,by=list(trial,cor,sub),mean))
names(trial_conf_sub) <- c("trial","cor","sub","cj")

pred_conf_sub_beta_learn <- with(Data_beta,aggregate(cj_pred_beta_learn,by=list(trial,cor,sub),mean))
names(pred_conf_sub_beta_learn) <- c("trial","cor","sub","cj")

pred_conf_sub_alpha_learn <- with(Data_beta,aggregate(cj_pred_alpha_learn,by=list(trial,cor,sub),mean))
names(pred_conf_sub_alpha_learn) <- c("trial","cor","sub","cj")

pred_conf_sub_both_learn <- with(Data_beta,aggregate(cj_pred_both_learn,by=list(trial,cor,sub),mean))
names(pred_conf_sub_both_learn) <- c("trial","cor","sub","cj")

pred_conf_sub_no_learn <- with(Data_beta,aggregate(cj_pred_no_learn,by=list(trial,cor,sub),mean))
names(pred_conf_sub_no_learn) <- c("trial","cor","sub","cj")

trials <- data.frame(trial=rep((0:(Ntrials-1)),each=2),
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
  if (s %in% subs_beta) {
    cj_pred_ma_beta_learn[cj_pred_ma_beta_learn$sub==s&cj_pred_ma_beta_learn$cor==0,"cj"] <- ma(subset(cj_pred_ma_beta_learn,sub==s&cor==0),n_err,"cj")
    cj_pred_ma_beta_learn[cj_pred_ma_beta_learn$sub==s&cj_pred_ma_beta_learn$cor==1,"cj"] <- ma(subset(cj_pred_ma_beta_learn,sub==s&cor==1),n,"cj")
    cj_pred_ma_alpha_learn[cj_pred_ma_alpha_learn$sub==s&cj_pred_ma_alpha_learn$cor==0,"cj"] <- ma(subset(cj_pred_ma_alpha_learn,sub==s&cor==0),n_err,"cj")
    cj_pred_ma_alpha_learn[cj_pred_ma_alpha_learn$sub==s&cj_pred_ma_alpha_learn$cor==1,"cj"] <- ma(subset(cj_pred_ma_alpha_learn,sub==s&cor==1),n,"cj")
    cj_pred_ma_both_learn[cj_pred_ma_both_learn$sub==s&cj_pred_ma_both_learn$cor==0,"cj"] <- ma(subset(cj_pred_ma_both_learn,sub==s&cor==0),n_err,"cj")
    cj_pred_ma_both_learn[cj_pred_ma_both_learn$sub==s&cj_pred_ma_both_learn$cor==1,"cj"] <- ma(subset(cj_pred_ma_both_learn,sub==s&cor==1),n,"cj")
    cj_pred_ma_no_learn[cj_pred_ma_no_learn$sub==s&cj_pred_ma_no_learn$cor==0,"cj"] <- ma(subset(cj_pred_ma_no_learn,sub==s&cor==0),n_err,"cj")
    cj_pred_ma_no_learn[cj_pred_ma_no_learn$sub==s&cj_pred_ma_no_learn$cor==1,"cj"] <- ma(subset(cj_pred_ma_no_learn,sub==s&cor==1),n,"cj")
  }
}

cj_pred_ma_beta_learn$model <- "beta"
cj_pred_ma_alpha_learn$model <- "alpha"
cj_pred_ma_no_learn$model <- "no"
cj_pred_ma_both_learn$model <- "both"
cj_pred <- rbind(cj_pred_ma_no_learn,cj_pred_ma_alpha_learn,
                 cj_pred_ma_beta_learn,cj_pred_ma_both_learn)
# Behavior results figure layout ----------------------------------------------------
diff_order <- c('hard','medium','easy')
Ndiff <- length(diff_order)
Nbeta <- length(unique(Data_beta$sub))
cexax <- 1
cex.lab <- 1
cexleg <- 1
cex.title <- 1
cex.axis <- 1
linelab <- 2.75
lwdmean <- 3
mar.rt.acc <- c(2,3.75,0,1)
mar.cj <- c(4,3.75,0,1)
mar.raw.cj <- c(4,3.75,1,1)
yrange_rt <- c(.8,1.1)
yrange_cj <- c(4.2,5.4)
col_plus_overlay <- rgb(213,94,0,51,maxColorValue = 255)
col_minus_overlay <- rgb(0,114,178,51,maxColorValue = 255)

setwd(paste0(curdir,"/plots"))
jpeg("Behavior_results2.jpg", width = 19, height = 24, units = "cm", pointsize = 15, res = 1000)

# layout(matrix(c(1,2,3,4,5,6),ncol=2,byrow = F),widths = c(1,2),heights = c(1,1,1.5))
layout(matrix(c(1,2,3,7,4,5,6,7),ncol=2,byrow = F),widths = c(1,2),heights = c(1,1,1.5,1.5))

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
     ylim=c(.5,1),xaxt='n',las=2)
title(ylab="Accuracy",xlab="", line = linelab, cex.lab = cex.lab)
axis(1,1:Ndiff,c("","",""),cex.axis=cexax)
points(colMeans(plot_cor_low),type='b',lwd=lwdmean,col=BLUE)
error.bar(1:Ndiff,colMeans(plot_cor_low),
          colSds(plot_cor_low,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=BLUE)

points(colMeans(plot_cor_high),type='b',lwd=lwdmean,col=VERMILLION)
error.bar(1:Ndiff,colMeans(plot_cor_high),
          colSds(plot_cor_high,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=VERMILLION)
legend("bottom",border=F,legend=c("Low","High"),pch=16,horiz=T,
       col=c(BLUE,VERMILLION),bty="n",
       cex=cexleg,title = "Feedback condition")

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
yrange <- c(.75,1)
plot(xlab="",ylab="",colMeans(plot_rt_low),frame=F,type='n',xlim=c(.8,Ndiff+.2),
     ylim=yrange,xaxt='n',las=2)
title(ylab="Reaction time (s)",xlab="", line = linelab, cex.lab = cex.lab)
# axis(1,1:Ndiff,c("Hard","Average","Easy"),cex.axis=cexax)
axis(1,1:Ndiff,c("","",""),cex.axis=cexax)
points(colMeans(plot_rt_low),type='b',lwd=lwdmean,col=BLUE)
error.bar(1:Ndiff,colMeans(plot_rt_low),
          colSds(plot_rt_low,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=BLUE)

points(colMeans(plot_rt_high),type='b',lwd=lwdmean,col=VERMILLION)
error.bar(1:Ndiff,colMeans(plot_rt_high),
          colSds(plot_rt_high,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=VERMILLION)

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
     ylim=yrange,xaxt='n',las=2)
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
legend(1.7,5,border=F,legend=c("Correct","Incorrect"),lty=c(1,2),horiz=F,
       col="black",bty="n",seg.len=1.5,
       cex=cexleg,title = "Trial accuracy")

# Plot aggregated traces  -------------------------------------------------
m <- 'both'

cj_pred <- with(subset(anal_sim,model==m),aggregate(cj_scaled,by=list(trial,sub),mean))
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

conf_group <- with(Data,aggregate(cj_integer,by=list(phase_block,manip,sub,condition,cor),mean))
names(conf_group) <- c('phase_block','manip','sub','condition','cor','cj')
conf_group <- merge(conf_group,trials_phase,all=T)
conf_group[conf_group$sub %in% Data_alpha$sub,'manip'] <- 'alpha'
conf_group[conf_group$sub %in% Data_beta$sub,'manip'] <- 'beta'
conf_beta <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_count <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,length)
conf_beta_sd <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_beta_sd[,2:5] <- conf_beta_sd[,2:5]/sqrt(conf_beta_count[,2:5])

xlen <- dim(conf_beta)[1]


par(mar=mar.rt.acc)

plot(cor_beta$minus,ylim=c(.5,1),col = BLUE, type = 'b',main="",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab,cex.axis=cex.axis,las=2)
axis(1, at = 1:Nphase_block, labels = F,cex.axis=cex.axis)
title(ylab = "Accuracy", xlab = "", 
      line = linelab,cex.lab=cex.lab)
lines(cor_beta$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
error.bar(1:xlen,cor_beta$minus,cor_beta_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:xlen,cor_beta$plus,cor_beta_sd$plus,
          lwd=2, col = VERMILLION)

par(mar=mar.rt.acc)

plot(rt_beta$minus,ylim=yrange_rt,col = BLUE, type = 'b',main="",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab,cex.axis=cex.axis,las=2)
axis(1, at = 1:Nphase_block, labels = F,cex.axis=cex.axis)
title(ylab = "Reaction time (s)", xlab = "", 
      line = linelab,cex.lab=cex.lab)
lines(rt_beta$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
error.bar(1:xlen,rt_beta$minus,rt_beta_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:xlen,rt_beta$plus,rt_beta_sd$plus,
          lwd=2, col = VERMILLION)

par(mar=mar.cj)


plot(conf_beta$minus_0,ylim=yrange_cj,col = BLUE, type = 'b',main="",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab,cex.axis=cex.axis,las=2)
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
        border=F,col=col_minus_overlay)
polygon(c(1:xlen,xlen:1),c(conf_beta_pred$minus_1 + conf_beta_pred_sd$minus_1,
                           (conf_beta_pred$minus_1 - conf_beta_pred_sd$minus_1)[xlen:1]),
        border=F,col=col_minus_overlay)
polygon(c(1:xlen,xlen:1),c(conf_beta_pred$plus_0 + conf_beta_pred_sd$plus_0,
                           (conf_beta_pred$plus_0 - conf_beta_pred_sd$plus_0)[xlen:1]),
        border=F,col=col_plus_overlay)
polygon(c(1:xlen,xlen:1),c(conf_beta_pred$plus_1 + conf_beta_pred_sd$plus_1,
                           (conf_beta_pred$plus_1 - conf_beta_pred_sd$plus_1)[xlen:1]),
        border=F,col=col_plus_overlay)
legend("top",legend = c('Behavior','Model'),lty = c(1,NA),col = c("black",rgb(0,0,0,.5)),
       border=NA,pch = c(16,15),horiz = T, bty = 'n',cex = cexleg,pt.cex=c(1,2))

# Plot confidence rolling mean over the course of the experiment ---------------------------------------------------
par(mar=mar.raw.cj)

Nphase <- 6

# Plot ExpB trace 
plus_first <- unique(subset(Data_beta,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_beta,phase==0&condition=='minus')$sub)
Ntrials_phase <- Ntrials/length(unique(Data$phase))

conf_min <- with(subset(cj_ma,sub %in% minus_first),aggregate(cj*6, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% minus_first),aggregate(cj*6, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% plus_first),aggregate(cj*6, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% plus_first),aggregate(cj*6, by=list(trial,cor),se,na.rm=T))
names(conf_plus_se) <- c("trial","cor","cj")

trials_phase1 <- c(seq(0,Ntrials_phase-1), 
                   seq(0,Ntrials_phase-1) + Ntrials_phase*2, 
                   seq(0,Ntrials_phase-1) + Ntrials_phase*4)

xlen <- dim(conf_plus)[1]/2

conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj
conf_plus_cor <- subset(conf_plus,cor==1)$cj
conf_plus_cor_se <- subset(conf_plus_se,cor==1)$cj


plot(conf_min_cor,bty='n',lty = 2,type='l',col="white",ylim=c(4.6,5.5),las=2,     
     main= NULL,cex.lab = cex.lab,cex.axis=cex.axis,xaxt='n',yaxt='n',xlab='',ylab='')
title(xlab = "Trial", ylab = "Confidence", line = linelab, cex.lab = cex.lab)
axis(1,at=seq(0,Ntrials,Ntrials/6),labels=seq(0,Ntrials,Ntrials/6),cex.axis=cex.axis)
axis(2,at=seq(4.6,5.5,.3),labels=seq(4.6,5.5,.3),cex.axis=cex.axis,las=2)

abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')

for (phase in 1:Nphase) {
  if (phase == 1) {
    trials_phase <- seq(Ntrials_phase*(phase-1)+1,Ntrials_phase*phase)
  } else {
    trials_phase <- seq(Ntrials_phase*(phase-1),Ntrials_phase*phase)
  }
  
  if (phase %% 2 == 1) {
    
    lines(trials_phase,conf_min_cor[trials_phase],col=BLUE,lty=2)
    lines(trials_phase,conf_plus_cor[trials_phase],col=VERMILLION)
    
    polygon(c(trials_phase,trials_phase[Ntrials_phase:1]),
            c((conf_min_cor + conf_min_cor_se)[trials_phase],
              (conf_min_cor - conf_min_cor_se)[trials_phase[Ntrials_phase:1]]),
            border=F,col=col_minus_overlay)
    polygon(c(trials_phase,trials_phase[Ntrials_phase:1]),
            c((conf_plus_cor + conf_plus_cor_se)[trials_phase],
              (conf_plus_cor - conf_plus_cor_se)[trials_phase[Ntrials_phase:1]]),
            border=F,col=col_plus_overlay)
  } else {
    
    lines(trials_phase,conf_min_cor[trials_phase],col=VERMILLION,lty=2)
    lines(trials_phase,conf_plus_cor[trials_phase],col=BLUE)
    
    polygon(c(trials_phase,trials_phase[Ntrials_phase:1]),
            c((conf_min_cor + conf_min_cor_se)[trials_phase],
              (conf_min_cor - conf_min_cor_se)[trials_phase[Ntrials_phase:1]]),
            border=F,col=col_plus_overlay)
    polygon(c(trials_phase,trials_phase[Ntrials_phase:1]),
            c((conf_plus_cor + conf_plus_cor_se)[trials_phase],
              (conf_plus_cor - conf_plus_cor_se)[trials_phase[Ntrials_phase:1]]),
            border=F,col=col_minus_overlay)
  }
}  

legend("bottom",border=F,legend=c("High first","Low first"),lty=c(1,2),horiz=T,
       col="black",bty="n",seg.len=1.5,
       cex=cexleg,title = "Counterbalance order")

dev.off()
par(mar=c(5,4,4,2)+.1)

# Reanalyze confidence according to group -----------------------------------------------------------
if (stat_tests) {
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
  
  # Same results if we look at the difference in BIC instead of groups
  cj.cond.acc.interaction.beta.bic <- lmer(cj ~ condition*cor*difflevel*withinphasetrial + condition:withinphasetrial:diff_learn + (cor + condition + condition:cor|sub),
                                           data = Data_beta, REML = F,control = control)
  anova(cj.cond.acc.interaction.beta.bic)
}


# Dynamic confidence plot according to group --------------------------------------------------------
trials_phase <- data.frame(phase_block=rep(0:(Nphase_block-1),each=2),
                           condition=c("minus","plus"),sub=rep(subs[subs%in%Data_beta$sub],each=Nphase_block*2))

conf_group <- with(Data_beta,aggregate(cj_integer,by=list(phase_block,sub,condition,cor,group),mean))
names(conf_group) <- c('phase_block','sub','condition','cor','group','cj')
conf_group <- merge(conf_group,trials_phase,all=T)
conf_static <- subset(conf_group,group=="static")
conf_dynamic <- subset(conf_group,group=="dynamic")

conf_beta_static <- cast(conf_static,phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_static_count <- cast(conf_static,phase_block~condition+cor,length)
conf_beta_static_sd <- cast(conf_static,phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_beta_static_sd[,2:5] <- conf_beta_static_sd[,2:5]/sqrt(conf_beta_static_count[,2:5])

conf_beta_dynamic <- cast(conf_dynamic,phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_dynamic_count <- cast(conf_dynamic,phase_block~condition+cor,length)
conf_beta_dynamic_sd <- cast(conf_dynamic,phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_beta_dynamic_sd[,2:5] <- conf_beta_dynamic_sd[,2:5]/sqrt(conf_beta_dynamic_count[,2:5])


cj_pred_static_nolearn <- with(subset(anal_sim,model=="no" & sub%in%sub_nolearn_best),aggregate(cj_scaled,by=list(trial,sub),mean))
names(cj_pred_static_nolearn) <- c("trial","sub","cj")
cj_pred_static_nolearn <- merge(cj_pred_static_nolearn,Data_beta[,c("trial","sub","phase_block","condition","cor")])

conf_group_pred <- with(cj_pred_static_nolearn,aggregate(cj,by=list(phase_block,sub,condition,cor),mean))
names(conf_group_pred) <- c('phase_block','sub','condition','cor','cj')

conf_beta_static_pred_nolearn_count <- cast(conf_group_pred,phase_block~condition+cor, fun.aggregate = length)
conf_beta_static_pred_nolearn <- cast(conf_group_pred,phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_static_pred_nolearn_sd <- cast(conf_group_pred,phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_beta_static_pred_nolearn_sd[,2:5] <- conf_beta_static_pred_nolearn_sd[,2:5]/sqrt(conf_beta_static_pred_nolearn_count[,2:5])

cj_pred_dynamic_nolearn <- with(subset(anal_sim,model=="no" & sub%in%sub_learn_best),aggregate(cj_scaled,by=list(trial,sub),mean))
names(cj_pred_dynamic_nolearn) <- c("trial","sub","cj")
cj_pred_dynamic_nolearn <- merge(cj_pred_dynamic_nolearn,Data_beta[,c("trial","sub","phase_block","condition","cor")])

conf_group_pred <- with(cj_pred_dynamic_nolearn,aggregate(cj,by=list(phase_block,sub,condition,cor),mean))
names(conf_group_pred) <- c('phase_block','sub','condition','cor','cj')

conf_beta_dynamic_pred_nolearn_count <- cast(conf_group_pred,phase_block~condition+cor, fun.aggregate = length)
conf_beta_dynamic_pred_nolearn <- cast(conf_group_pred,phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_dynamic_pred_nolearn_sd <- cast(conf_group_pred,phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_beta_dynamic_pred_nolearn_sd[,2:5] <- conf_beta_dynamic_pred_nolearn_sd[,2:5]/sqrt(conf_beta_dynamic_pred_nolearn_count[,2:5])


cj_pred_static_learn <- with(subset(anal_sim,model=="both" & sub%in%sub_nolearn_best),aggregate(cj_scaled,by=list(trial,sub),mean))
names(cj_pred_static_learn) <- c("trial","sub","cj")
cj_pred_static_learn <- merge(cj_pred_static_learn,Data_beta[,c("trial","sub","phase_block","condition","cor")])

conf_group_pred <- with(cj_pred_static_learn,aggregate(cj,by=list(phase_block,sub,condition,cor),mean))
names(conf_group_pred) <- c('phase_block','sub','condition','cor','cj')

conf_beta_static_pred_learn_count <- cast(conf_group_pred,phase_block~condition+cor, fun.aggregate = length)
conf_beta_static_pred_learn <- cast(conf_group_pred,phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_static_pred_learn_sd <- cast(conf_group_pred,phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_beta_static_pred_learn_sd[,2:5] <- conf_beta_static_pred_learn_sd[,2:5]/sqrt(conf_beta_static_pred_learn_count[,2:5])


cj_pred_dynamic_learn <- with(subset(anal_sim,model=="both" & sub%in%sub_learn_best),aggregate(cj_scaled,by=list(trial,sub),mean))
names(cj_pred_dynamic_learn) <- c("trial","sub","cj")
cj_pred_dynamic_learn <- merge(cj_pred_dynamic_learn,Data_beta[,c("trial","sub","phase_block","condition","cor")])

conf_group_pred <- with(cj_pred_dynamic_learn,aggregate(cj,by=list(phase_block,sub,condition,cor),mean))
names(conf_group_pred) <- c('phase_block','sub','condition','cor','cj')

conf_beta_dynamic_pred_learn_count <- cast(conf_group_pred,phase_block~condition+cor, fun.aggregate = length)
conf_beta_dynamic_pred_learn <- cast(conf_group_pred,phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_dynamic_pred_learn_sd <- cast(conf_group_pred,phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_beta_dynamic_pred_learn_sd[,2:5] <- conf_beta_dynamic_pred_learn_sd[,2:5]/sqrt(conf_beta_dynamic_pred_learn_count[,2:5])

cex_size <- function(size,cex.layout) {
  return(size/(par()$ps*cex.layout))
}
### Adjust sizes and positions
cex.layout <- .83
# Text
cex.legend <- cex_size(8,cex.layout)
cex.axis <- cex_size(8,cex.layout) 
cex.lab <- cex_size(10,cex.layout)
line.width <- 1
cex.pt <- .8
yrange <- c(4,5.5)

xlen <- dim(conf_beta_static)[1]

jpeg("conf_beta_group_new.jpeg",width=16.5,height=11,units="cm",res=600)
par(mfrow=c(2,2))
par(mar=c(3,3,0,0)+.1)

# dynamic group, learn model
plot(conf_beta_dynamic$minus_0,ylim=yrange,col = BLUE, type = 'l',main="",
     lty = 2,  lwd = line.width, bty = 'n', xaxt = 'n', xlab="",ylab="",cex.axis=cex.axis)
title(ylab = "Confidence",xlab = "",cex.lab=cex.lab,line=2)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis)
lines(conf_beta_dynamic$plus_0, type = 'l',  col = VERMILLION, lwd = line.width, lty = 2)
lines(conf_beta_dynamic$plus_1, type = 'l',  col = VERMILLION, lwd = line.width, lty = 1)
lines(conf_beta_dynamic$minus_1, type = 'l',  col = BLUE, lwd = line.width, lty = 1)
points(conf_beta_dynamic$plus_0, pch = 16, cex = cex.pt, col = VERMILLION)
points(conf_beta_dynamic$plus_1, pch = 16, cex = cex.pt, col = VERMILLION)
points(conf_beta_dynamic$minus_1, pch = 16, cex = cex.pt, col = BLUE)
points(conf_beta_dynamic$minus_0, pch = 16, cex = cex.pt, col = BLUE)
error.bar(1:length(conf_beta_dynamic$minus_0),conf_beta_dynamic$minus_0,conf_beta_dynamic_sd$minus_0,
          lwd = line.width, col = BLUE)
error.bar(1:length(conf_beta_dynamic$plus_0),conf_beta_dynamic$plus_0,conf_beta_dynamic_sd$plus_0,
          lwd = line.width, col = VERMILLION)
error.bar(1:length(conf_beta_dynamic$minus_1),conf_beta_dynamic$minus_1,conf_beta_dynamic_sd$minus_1,
          lwd = line.width, col = BLUE)
error.bar(1:length(conf_beta_dynamic$plus_1),conf_beta_dynamic$plus_1,conf_beta_dynamic_sd$plus_1,
          lwd = line.width, col = VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_beta_dynamic_pred_learn$minus_0 + conf_beta_dynamic_pred_learn_sd$minus_0,
                           (conf_beta_dynamic_pred_learn$minus_0 - conf_beta_dynamic_pred_learn_sd$minus_0)[xlen:1]),
        border=F,col=col_minus_overlay)
polygon(c(1:xlen,xlen:1),c(conf_beta_dynamic_pred_learn$minus_1 + conf_beta_dynamic_pred_learn_sd$minus_1,
                           (conf_beta_dynamic_pred_learn$minus_1 - conf_beta_dynamic_pred_learn_sd$minus_1)[xlen:1]),
        border=F,col=col_minus_overlay)
polygon(c(1:xlen,xlen:1),c(conf_beta_dynamic_pred_learn$plus_0 + conf_beta_dynamic_pred_learn_sd$plus_0,
                           (conf_beta_dynamic_pred_learn$plus_0 - conf_beta_dynamic_pred_learn_sd$plus_0)[xlen:1]),
        border=F,col=col_plus_overlay)
polygon(c(1:xlen,xlen:1),c(conf_beta_dynamic_pred_learn$plus_1 + conf_beta_dynamic_pred_learn_sd$plus_1,
                           (conf_beta_dynamic_pred_learn$plus_1 - conf_beta_dynamic_pred_learn_sd$plus_1)[xlen:1]),
        border=F,col=col_plus_overlay)
legend(0.5,.93*6,legend = c('High feedback','Low feedback'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16), bty = 'n',cex = cex.legend)
legend(11.5,.93*6,legend=c("Correct trials","Error trials"),
       title = NULL,lty=c(1,2),bty = "n",inset=0,
       cex = cex.legend, lwd=lwd.dat,seg.len=1.5)

# dynamic group, no learn model
plot(conf_beta_dynamic$minus_0,ylim=yrange,col = BLUE, type = 'l',main="",
     lty = 2,  lwd = line.width, bty = 'n', xaxt = 'n', xlab="",ylab="",cex.axis=cex.axis)
title(ylab = "",xlab = "",cex.lab=cex.lab,line=2)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis)
lines(conf_beta_dynamic$plus_0, type = 'l',  col = VERMILLION, lwd = line.width, lty = 2)
lines(conf_beta_dynamic$plus_1, type = 'l',  col = VERMILLION, lwd = line.width, lty = 1)
lines(conf_beta_dynamic$minus_1, type = 'l',  col = BLUE, lwd = line.width, lty = 1)
points(conf_beta_dynamic$plus_0, pch = 16, cex = cex.pt, col = VERMILLION)
points(conf_beta_dynamic$plus_1, pch = 16, cex = cex.pt, col = VERMILLION)
points(conf_beta_dynamic$minus_1, pch = 16, cex = cex.pt, col = BLUE)
points(conf_beta_dynamic$minus_0, pch = 16, cex = cex.pt, col = BLUE)
error.bar(1:length(conf_beta_dynamic$minus_0),conf_beta_dynamic$minus_0,conf_beta_dynamic_sd$minus_0,
          lwd = line.width, col = BLUE)
error.bar(1:length(conf_beta_dynamic$plus_0),conf_beta_dynamic$plus_0,conf_beta_dynamic_sd$plus_0,
          lwd = line.width, col = VERMILLION)
error.bar(1:length(conf_beta_dynamic$minus_1),conf_beta_dynamic$minus_1,conf_beta_dynamic_sd$minus_1,
          lwd = line.width, col = BLUE)
error.bar(1:length(conf_beta_dynamic$plus_1),conf_beta_dynamic$plus_1,conf_beta_dynamic_sd$plus_1,
          lwd = line.width, col = VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_beta_dynamic_pred_nolearn$minus_0 + conf_beta_dynamic_pred_nolearn_sd$minus_0,
                           (conf_beta_dynamic_pred_nolearn$minus_0 - conf_beta_dynamic_pred_nolearn_sd$minus_0)[xlen:1]),
        border=F,col=col_minus_overlay)
polygon(c(1:xlen,xlen:1),c(conf_beta_dynamic_pred_nolearn$minus_1 + conf_beta_dynamic_pred_nolearn_sd$minus_1,
                           (conf_beta_dynamic_pred_nolearn$minus_1 - conf_beta_dynamic_pred_nolearn_sd$minus_1)[xlen:1]),
        border=F,col=col_minus_overlay)
polygon(c(1:xlen,xlen:1),c(conf_beta_dynamic_pred_nolearn$plus_0 + conf_beta_dynamic_pred_nolearn_sd$plus_0,
                           (conf_beta_dynamic_pred_nolearn$plus_0 - conf_beta_dynamic_pred_nolearn_sd$plus_0)[xlen:1]),
        border=F,col=col_plus_overlay)
polygon(c(1:xlen,xlen:1),c(conf_beta_dynamic_pred_nolearn$plus_1 + conf_beta_dynamic_pred_nolearn_sd$plus_1,
                           (conf_beta_dynamic_pred_nolearn$plus_1 - conf_beta_dynamic_pred_nolearn_sd$plus_1)[xlen:1]),
        border=F,col=col_plus_overlay)
legend("top",legend = c('Behavior','Model'),lty = c(1,NA),col = c("black",rgb(0,0,0,.5)),
       border=NA,pch = c(16,15),horiz = T, bty = 'n',cex = cex.legend,pt.cex=c(1,2))

# Static group, learn model
plot(conf_beta_static$minus_0,ylim=yrange,col = BLUE, type = 'l',main="",
     lty = 2,  lwd = line.width, bty = 'n', xaxt = 'n', xlab="",ylab="",cex.axis=cex.axis)
title(ylab = "Confidence",xlab = "Within phase groups of 7 trials",cex.lab=cex.lab,line=2)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis)
lines(conf_beta_static$plus_0, type = 'l',  col = VERMILLION, lwd = line.width, lty = 2)
lines(conf_beta_static$plus_1, type = 'l',  col = VERMILLION, lwd = line.width, lty = 1)
lines(conf_beta_static$minus_1, type = 'l',  col = BLUE, lwd = line.width, lty = 1)
points(conf_beta_static$plus_0, pch = 16, cex = cex.pt, col = VERMILLION)
points(conf_beta_static$plus_1, pch = 16, cex = cex.pt, col = VERMILLION)
points(conf_beta_static$minus_1, pch = 16, cex = cex.pt, col = BLUE)
points(conf_beta_static$minus_0, pch = 16, cex = cex.pt, col = BLUE)
error.bar(1:length(conf_beta_static$minus_0),conf_beta_static$minus_0,conf_beta_static_sd$minus_0,
          lwd = line.width, col = BLUE)
error.bar(1:length(conf_beta_static$plus_0),conf_beta_static$plus_0,conf_beta_static_sd$plus_0,
          lwd = line.width, col = VERMILLION)
error.bar(1:length(conf_beta_static$minus_1),conf_beta_static$minus_1,conf_beta_static_sd$minus_1,
          lwd = line.width, col = BLUE)
error.bar(1:length(conf_beta_static$plus_1),conf_beta_static$plus_1,conf_beta_static_sd$plus_1,
          lwd = line.width, col = VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_beta_static_pred_learn$minus_0 + conf_beta_static_pred_learn_sd$minus_0,
                           (conf_beta_static_pred_learn$minus_0 - conf_beta_static_pred_learn_sd$minus_0)[xlen:1]),
        border=F,col=col_minus_overlay)
polygon(c(1:xlen,xlen:1),c(conf_beta_static_pred_learn$minus_1 + conf_beta_static_pred_learn_sd$minus_1,
                           (conf_beta_static_pred_learn$minus_1 - conf_beta_static_pred_learn_sd$minus_1)[xlen:1]),
        border=F,col=col_minus_overlay)
polygon(c(1:xlen,xlen:1),c(conf_beta_static_pred_learn$plus_0 + conf_beta_static_pred_learn_sd$plus_0,
                           (conf_beta_static_pred_learn$plus_0 - conf_beta_static_pred_learn_sd$plus_0)[xlen:1]),
        border=F,col=col_plus_overlay)
polygon(c(1:xlen,xlen:1),c(conf_beta_static_pred_learn$plus_1 + conf_beta_static_pred_learn_sd$plus_1,
                           (conf_beta_static_pred_learn$plus_1 - conf_beta_static_pred_learn_sd$plus_1)[xlen:1]),
        border=F,col=col_plus_overlay)

# Static group, no learn model
plot(conf_beta_static$minus_0,ylim=yrange,col = BLUE, type = 'l',main="",
     lty = 2,  lwd = line.width, bty = 'n', xaxt = 'n', xlab="",ylab="",cex.axis=cex.axis)
title(ylab = "",xlab = "Within phase groups of 7 trials",cex.lab=cex.lab,line=2)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis)
lines(conf_beta_static$plus_0, type = 'l',  col = VERMILLION, lwd = line.width, lty = 2)
lines(conf_beta_static$plus_1, type = 'l',  col = VERMILLION, lwd = line.width, lty = 1)
lines(conf_beta_static$minus_1, type = 'l',  col = BLUE, lwd = line.width, lty = 1)
points(conf_beta_static$plus_0, pch = 16, cex = cex.pt, col = VERMILLION)
points(conf_beta_static$plus_1, pch = 16, cex = cex.pt, col = VERMILLION)
points(conf_beta_static$minus_1, pch = 16, cex = cex.pt, col = BLUE)
points(conf_beta_static$minus_0, pch = 16, cex = cex.pt, col = BLUE)
error.bar(1:length(conf_beta_static$minus_0),conf_beta_static$minus_0,conf_beta_static_sd$minus_0,
          lwd = line.width, col = BLUE)
error.bar(1:length(conf_beta_static$plus_0),conf_beta_static$plus_0,conf_beta_static_sd$plus_0,
          lwd = line.width, col = VERMILLION)
error.bar(1:length(conf_beta_static$minus_1),conf_beta_static$minus_1,conf_beta_static_sd$minus_1,
          lwd = line.width, col = BLUE)
error.bar(1:length(conf_beta_static$plus_1),conf_beta_static$plus_1,conf_beta_static_sd$plus_1,
          lwd = line.width, col = VERMILLION)
polygon(c(1:xlen,xlen:1),c(conf_beta_static_pred_nolearn$minus_0 + conf_beta_static_pred_nolearn_sd$minus_0,
                           (conf_beta_static_pred_nolearn$minus_0 - conf_beta_static_pred_nolearn_sd$minus_0)[xlen:1]),
        border=F,col=col_minus_overlay)
polygon(c(1:xlen,xlen:1),c(conf_beta_static_pred_nolearn$minus_1 + conf_beta_static_pred_nolearn_sd$minus_1,
                           (conf_beta_static_pred_nolearn$minus_1 - conf_beta_static_pred_nolearn_sd$minus_1)[xlen:1]),
        border=F,col=col_minus_overlay)
polygon(c(1:xlen,xlen:1),c(conf_beta_static_pred_nolearn$plus_0 + conf_beta_static_pred_nolearn_sd$plus_0,
                           (conf_beta_static_pred_nolearn$plus_0 - conf_beta_static_pred_nolearn_sd$plus_0)[xlen:1]),
        border=F,col=col_plus_overlay)
polygon(c(1:xlen,xlen:1),c(conf_beta_static_pred_nolearn$plus_1 + conf_beta_static_pred_nolearn_sd$plus_1,
                           (conf_beta_static_pred_nolearn$plus_1 - conf_beta_static_pred_nolearn_sd$plus_1)[xlen:1]),
        border=F,col=col_plus_overlay)

dev.off()


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
# Better way to replace values because does not get an error if value to replace is not found
questions$questionResp <- ifelse(questions$questionResp=="['j']",5,questions$questionResp) 
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

ordinal.reg <- clm(data=subset(question2,sub %in% Data_beta$sub),questionResp~group)
anova(ordinal.reg)


# 3. 'Indicate to what extent you feel that the feedback you received closely 
# tracked your actual performance.'
question3 <- subset(questions,questionID==3)
question3$questionResp <- as.numeric(question3$questionResp)
with(question3,aggregate(questionResp,by=list(group,manip),table))
# Test difference in response between groups
question3$questionResp <- as.ordered(question3$questionResp)

ordinal.reg <- clm(data=subset(question3,sub %in% Data_beta$sub),questionResp~group)
anova(ordinal.reg)

# 4. Check participants understood the feedback
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

with(Data_beta,aggregate(sub,by=list(fb_understood,group),length_unique))

fb_understood.test.beta <- glmer(fb_understood ~ group + (1|sub), data=Data_beta, family="binomial", control = glmercontrol)
Anova(fb_understood.test.beta)

Data_alpha$fb_understood <- 0
Data_alpha[Data_alpha$sub %in% fb_understood_sub,]$fb_understood <- 1

with(Data_alpha,aggregate(sub,by=list(fb_understood,group),length_unique))

fb_understood.test.alpha <- glmer(fb_understood ~ group + (1|sub), data=Data_alpha, family="binomial", control = glmercontrol)
Anova(fb_understood.test.alpha)


# Plot Feedback -----------------------------------------------------------

diff_order <- c('hard','medium','easy')
Ndiff <- length(diff_order)
Nbeta <- length(unique(Data_beta$sub))
ptsize <- 10
cexax <- 1
cex.lab <- 1
cexleg <- 1*8/ptsize
cex.title <- 1
cex.axis <- 1
linelab <- 2.5
lwdmean <- 2
pch <- 16
mar.fb <- c(3.5,3.5,0,1)
jpeg("feedback.jpg", width = 6, height = 6, units = "cm", pointsize = ptsize, res = 1000)

par(mar=mar.fb)

plot_fb_low <- with(subset(Data_beta,condition=="minus"),
                    aggregate(fb*100,by=list(sub=sub,difflevel=difflevel,cor=cor),mean))
plot_fb_low_cor <- cast(subset(plot_fb_low,cor==1),sub~difflevel)
plot_fb_low_err <- cast(subset(plot_fb_low,cor==0),sub~difflevel)
plot_fb_low_cor <- plot_fb_low_cor[,diff_order] #Reorder columns to have hard -> easy
plot_fb_low_err <- plot_fb_low_err[,diff_order] #Reorder columns to have hard -> easy

plot_fb_high <- with(subset(Data_beta,condition=="plus"),
                     aggregate(fb*100,by=list(sub=sub,difflevel=difflevel,cor=cor),mean))
plot_fb_high_cor <- cast(subset(plot_fb_high,cor==1),sub~difflevel)
plot_fb_high_err <- cast(subset(plot_fb_high,cor==0),sub~difflevel)
plot_fb_high_cor <- plot_fb_high_cor[,diff_order] #Reorder columns to have hard -> easy
plot_fb_high_err <- plot_fb_high_err[,diff_order] #Reorder columns to have hard -> easy

yrange <- c(min(colMeans(rbind(plot_fb_low_err,plot_fb_high_err),na.rm = T)) - min(colSds(rbind(plot_fb_low_err,plot_fb_high_err),na.rm = T))/sqrt(Nbeta),
            max(colMeans(rbind(plot_fb_low_cor,plot_fb_high_cor),na.rm = T)) + max(colSds(rbind(plot_fb_low_cor,plot_fb_high_cor),na.rm = T))/sqrt(Nbeta))
plot(xlab="",ylab="",colMeans(plot_fb_low_cor),frame=F,type='n',xlim=c(.8,Ndiff+.2),
     ylim=c(0,100),xaxt='n',las=2)
title(ylab="Feedback",xlab="Trial difficulty", line = linelab, cex.lab = cex.lab)
axis(1,1:Ndiff,c("Hard","Average","Easy"),cex.axis=cexax)
# Incorrect trials
points(colMeans(plot_fb_low_err,na.rm=T),type='b',lwd=lwdmean,col=BLUE,lty=2,pch=pch)
error.bar(1:Ndiff,colMeans(plot_fb_low_err),
          colSds(plot_fb_low_err,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=BLUE,lty=1)
points(colMeans(plot_fb_high_err,na.rm=T),type='b',lwd=lwdmean,col=VERMILLION,lty=2,pch=pch)
error.bar(1:Ndiff,colMeans(plot_fb_high_err,na.rm=T),
          colSds(plot_fb_high_err,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=VERMILLION,lty=1)
# Correct trials
points(colMeans(plot_fb_low_cor,na.rm=T),type='b',lwd=lwdmean,col=BLUE,lty=1,pch=pch)
error.bar(1:Ndiff,colMeans(plot_fb_low_cor),
          colSds(plot_fb_low_cor,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=BLUE,lty=1)
points(colMeans(plot_fb_high_cor,na.rm=T),type='b',lwd=lwdmean,col=VERMILLION,lty=1,pch=pch)
error.bar(1:Ndiff,colMeans(plot_fb_high_cor,na.rm=T),
          colSds(plot_fb_high_cor,na.rm=T)/sqrt(Nbeta),lwd=lwdmean,length=0,col=VERMILLION,lty=1)
legend("bottomleft",border=F,legend=c("Correct","Incorrect"),lty=c(1,2),horiz=F,
       col="black",bty="n",seg.len=1.5,
       cex=cexleg,title = "Trial accuracy")
legend("bottomright",border=F,legend=c("High","Low"),lty=1,horiz=F,
       col=c(VERMILLION,BLUE),bty="n",seg.len=1.5,lwd=lwdmean,
       cex=cexleg,title = "Condition")

dev.off()

# Plot behavior results for alpha experiment ------------------------------
# Behavior results figure layout ----------------------------------------------------
diff_order <- c('hard','medium','easy')
Ndiff <- length(diff_order)
Nalpha <- length(unique(Data_alpha$sub))
cexax <- 1
cex.lab <- 1
cexleg <- 1
cex.title <- 1
cex.axis <- 1
linelab <- 2.75
lwdmean <- 3
mar.rt.acc <- c(2,3.75,0,1)
mar.cj <- c(4,3.75,0,1)
mar.raw.cj <- c(4,3.75,1,1)
yrange_rt <- c(.8,1.1)
yrange_cj <- c(4.2,5.4)
col_plus_overlay <- rgb(213,94,0,51,maxColorValue = 255)
col_minus_overlay <- rgb(0,114,178,51,maxColorValue = 255)

setwd(paste0(curdir,"/plots"))
jpeg("Behavior_results_alpha.jpg", width = 19, height = 24, units = "cm", pointsize = 15, res = 1000)

# layout(matrix(c(1,2,3,4,5,6),ncol=2,byrow = F),widths = c(1,2),heights = c(1,1,1.5))
layout(matrix(c(1,2,3,7,4,5,6,7),ncol=2,byrow = F),widths = c(1,2),heights = c(1,1,1.5,1.5))

###
# Accuracy
###
par(mar=mar.rt.acc)


plot_cor_low <- with(subset(Data_alpha,condition=="minus"),
                     aggregate(cor,by=list(sub=sub,difflevel=difflevel),mean))
plot_cor_low <- cast(plot_cor_low,sub~difflevel)
plot_cor_low <- plot_cor_low[,diff_order] #Reorder columns to have hard -> easy

plot_cor_high <- with(subset(Data_alpha,condition=="plus"),
                      aggregate(cor,by=list(sub=sub,difflevel=difflevel),mean))
plot_cor_high <- cast(plot_cor_high,sub~difflevel)
plot_cor_high <- plot_cor_high[,diff_order] #Reorder columns to have hard -> easy

yrange <- c(min(colMeans(rbind(plot_cor_low,plot_cor_high))) - min(colSds(rbind(plot_cor_low,plot_cor_high)))/sqrt(Nalpha),
            max(colMeans(rbind(plot_cor_low,plot_cor_high))) + max(colSds(rbind(plot_cor_low,plot_cor_high)))/sqrt(Nalpha))
plot(xlab="",ylab="",colMeans(plot_cor_low),frame=F,type='n',xlim=c(.8,Ndiff+.2),
     ylim=c(.5,1),xaxt='n',las=2)
title(ylab="Accuracy",xlab="", line = linelab, cex.lab = cex.lab)
axis(1,1:Ndiff,c("","",""),cex.axis=cexax)
points(colMeans(plot_cor_low),type='b',lwd=lwdmean,col=BLUE)
error.bar(1:Ndiff,colMeans(plot_cor_low),
          colSds(plot_cor_low,na.rm=T)/sqrt(Nalpha),lwd=lwdmean,length=0,col=BLUE)

points(colMeans(plot_cor_high),type='b',lwd=lwdmean,col=VERMILLION)
error.bar(1:Ndiff,colMeans(plot_cor_high),
          colSds(plot_cor_high,na.rm=T)/sqrt(Nalpha),lwd=lwdmean,length=0,col=VERMILLION)
legend("bottom",border=F,legend=c("Low","High"),pch=16,horiz=T,
       col=c(BLUE,VERMILLION),bty="n",
       cex=cexleg,title = "Feedback condition")

par(mar=mar.rt.acc)

###
# RT
###

plot_rt_low <- with(subset(Data_alpha,condition=="minus"&cor==1),
                    aggregate(rt,by=list(sub=sub,difflevel=difflevel),mean))
plot_rt_low <- cast(plot_rt_low,sub~difflevel)
plot_rt_low <- plot_rt_low[,diff_order] #Reorder columns to have hard -> easy

plot_rt_high <- with(subset(Data_alpha,condition=="plus"&cor==1),
                     aggregate(rt,by=list(sub=sub,difflevel=difflevel),mean))
plot_rt_high <- cast(plot_rt_high,sub~difflevel)
plot_rt_high <- plot_rt_high[,diff_order] #Reorder columns to have hard -> easy

yrange <- c(min(c(colMeans(plot_rt_low),colMeans(plot_rt_high))) - max(c(colSds(plot_rt_low),colSds(plot_rt_high)))/sqrt(Nalpha),
            max(c(colMeans(plot_rt_low),colMeans(plot_rt_high))) + max(c(colSds(plot_rt_low),colSds(plot_rt_high)))/sqrt(Nalpha))
yrange <- c(.75,1)
plot(xlab="",ylab="",colMeans(plot_rt_low),frame=F,type='n',xlim=c(.8,Ndiff+.2),
     ylim=yrange,xaxt='n',las=2)
title(ylab="Reaction time (s)",xlab="", line = linelab, cex.lab = cex.lab)
# axis(1,1:Ndiff,c("Hard","Average","Easy"),cex.axis=cexax)
axis(1,1:Ndiff,c("","",""),cex.axis=cexax)
points(colMeans(plot_rt_low),type='b',lwd=lwdmean,col=BLUE)
error.bar(1:Ndiff,colMeans(plot_rt_low),
          colSds(plot_rt_low,na.rm=T)/sqrt(Nalpha),lwd=lwdmean,length=0,col=BLUE)

points(colMeans(plot_rt_high),type='b',lwd=lwdmean,col=VERMILLION)
error.bar(1:Ndiff,colMeans(plot_rt_high),
          colSds(plot_rt_high,na.rm=T)/sqrt(Nalpha),lwd=lwdmean,length=0,col=VERMILLION)

###
# Confidence
###
par(mar=mar.cj)

plot_cj_low <- with(subset(Data_alpha,condition=="minus"),
                    aggregate(cj_integer,by=list(sub=sub,difflevel=difflevel,cor=cor),mean))
plot_cj_low_cor <- cast(subset(plot_cj_low,cor==1),sub~difflevel)
plot_cj_low_err <- cast(subset(plot_cj_low,cor==0),sub~difflevel)
plot_cj_low_cor <- plot_cj_low_cor[,diff_order] #Reorder columns to have hard -> easy
plot_cj_low_err <- plot_cj_low_err[,diff_order] #Reorder columns to have hard -> easy

plot_cj_high <- with(subset(Data_alpha,condition=="plus"),
                     aggregate(cj_integer,by=list(sub=sub,difflevel=difflevel,cor=cor),mean))
plot_cj_high_cor <- cast(subset(plot_cj_high,cor==1),sub~difflevel)
plot_cj_high_err <- cast(subset(plot_cj_high,cor==0),sub~difflevel)
plot_cj_high_cor <- plot_cj_high_cor[,diff_order] #Reorder columns to have hard -> easy
plot_cj_high_err <- plot_cj_high_err[,diff_order] #Reorder columns to have hard -> easy

yrange <- c(min(colMeans(rbind(plot_cj_low_err,plot_cj_high_err),na.rm = T)) - min(colSds(rbind(plot_cj_low_err,plot_cj_high_err),na.rm = T))/sqrt(Nalpha),
            max(colMeans(rbind(plot_cj_low_cor,plot_cj_high_cor),na.rm = T)) + max(colSds(rbind(plot_cj_low_cor,plot_cj_high_cor),na.rm = T))/sqrt(Nalpha))
plot(xlab="",ylab="",colMeans(plot_cj_low_cor),frame=F,type='n',xlim=c(.8,Ndiff+.2),
     ylim=yrange,xaxt='n',las=2)
title(ylab="Confidence",xlab="Trial difficulty", line = linelab, cex.lab = cex.lab)
axis(1,1:Ndiff,c("Hard","Average","Easy"),cex.axis=cexax)
# Incorrect trials
points(colMeans(plot_cj_low_err,na.rm=T),type='b',lwd=lwdmean,col=BLUE,lty=2)
error.bar(1:Ndiff,colMeans(plot_cj_low_err),
          colSds(plot_cj_low_err,na.rm=T)/sqrt(Nalpha),lwd=lwdmean,length=0,col=BLUE,lty=1)
points(colMeans(plot_cj_high_err,na.rm=T),type='b',lwd=lwdmean,col=VERMILLION,lty=2)
error.bar(1:Ndiff,colMeans(plot_cj_high_err,na.rm=T),
          colSds(plot_cj_high_err,na.rm=T)/sqrt(Nalpha),lwd=lwdmean,length=0,col=VERMILLION,lty=1)
# Correct trials
points(colMeans(plot_cj_low_cor,na.rm=T),type='b',lwd=lwdmean,col=BLUE,lty=1)
error.bar(1:Ndiff,colMeans(plot_cj_low_cor),
          colSds(plot_cj_low_cor,na.rm=T)/sqrt(Nalpha),lwd=lwdmean,length=0,col=BLUE,lty=1)
points(colMeans(plot_cj_high_cor,na.rm=T),type='b',lwd=lwdmean,col=VERMILLION,lty=1)
error.bar(1:Ndiff,colMeans(plot_cj_high_cor,na.rm=T),
          colSds(plot_cj_high_cor,na.rm=T)/sqrt(Nalpha),lwd=lwdmean,length=0,col=VERMILLION,lty=1)
legend(1.7,5,border=F,legend=c("Correct","Incorrect"),lty=c(1,2),horiz=F,
       col="black",bty="n",seg.len=1.5,
       cex=cexleg,title = "Trial accuracy")

# Plot aggregated traces  -------------------------------------------------
if (is.factor(Data$cor)) {
  Data$cor <- as.numeric(Data$cor)-1
}


trials_phase <- data.frame(phase_block=rep(0:(Nphase_block-1),each=2),
                           condition=c("minus","plus"),sub=rep(subs,each=Nphase_block*2))

rt_group <- with(Data,aggregate(rt,by=list(phase_block,manip,sub,condition),mean))
names(rt_group) <- c('phase_block','manip','sub','condition','cj')
rt_group <- merge(rt_group,trials_phase,all=T)
rt_group[rt_group$sub %in% Data_alpha$sub,'manip'] <- 'alpha'

rt_alpha <- cast(subset(rt_group,manip=='alpha'),phase_block~condition,fun.aggregate = mean,na.rm=T)
rt_alpha_count <- cast(subset(rt_group,manip=='alpha'),phase_block~condition,fun.aggregate=length)
rt_alpha_sd <- cast(subset(rt_group,manip=='alpha'),phase_block~condition,fun.aggregate = sd,na.rm=T)
rt_alpha_sd[,2:3] <- rt_alpha_sd[,2:3]/sqrt(rt_alpha_count[,2:3])

cor_group <- with(Data,aggregate(cor,by=list(phase_block,manip,sub,condition),mean))
names(cor_group) <- c('phase_block','manip','sub','condition','cj')
cor_group <- merge(cor_group,trials_phase,all=T)
cor_group[cor_group$sub %in% Data_alpha$sub,'manip'] <- 'alpha'


cor_alpha <- cast(subset(cor_group,manip=='alpha'),phase_block~condition,fun.aggregate = mean,na.rm=T)
cor_alpha_count <- cast(subset(cor_group,manip=='alpha'),phase_block~condition,fun.aggregate=length)
cor_alpha_sd <- cast(subset(cor_group,manip=='alpha'),phase_block~condition,fun.aggregate = sd,na.rm=T)
cor_alpha_sd[,2:3] <- cor_alpha_sd[,2:3]/sqrt(cor_alpha_count[,2:3])

conf_group <- with(Data,aggregate(cj_integer,by=list(phase_block,manip,sub,condition,cor),mean))
names(conf_group) <- c('phase_block','manip','sub','condition','cor','cj')
conf_group <- merge(conf_group,trials_phase,all=T)
conf_group[conf_group$sub %in% Data_alpha$sub,'manip'] <- 'alpha'
conf_group[conf_group$sub %in% Data_alpha$sub,'manip'] <- 'alpha'
conf_alpha <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_alpha_count <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor,length)
conf_alpha_sd <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor,fun.aggregate = sd,na.rm=T)
conf_alpha_sd[,2:5] <- conf_alpha_sd[,2:5]/sqrt(conf_alpha_count[,2:5])

xlen <- dim(conf_alpha)[1]


par(mar=mar.rt.acc)

plot(cor_alpha$minus,ylim=c(.5,1),col = BLUE, type = 'b',main="",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab,cex.axis=cex.axis,las=2)
axis(1, at = 1:Nphase_block, labels = F,cex.axis=cex.axis)
title(ylab = "Accuracy", xlab = "", 
      line = linelab,cex.lab=cex.lab)
lines(cor_alpha$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
error.bar(1:xlen,cor_alpha$minus,cor_alpha_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:xlen,cor_alpha$plus,cor_alpha_sd$plus,
          lwd=2, col = VERMILLION)

par(mar=mar.rt.acc)

plot(rt_alpha$minus,ylim=yrange_rt,col = BLUE, type = 'b',main="",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab,cex.axis=cex.axis,las=2)
axis(1, at = 1:Nphase_block, labels = F,cex.axis=cex.axis)
title(ylab = "Reaction time (s)", xlab = "", 
      line = linelab,cex.lab=cex.lab)
lines(rt_alpha$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
error.bar(1:xlen,rt_alpha$minus,rt_alpha_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:xlen,rt_alpha$plus,rt_alpha_sd$plus,
          lwd=2, col = VERMILLION)

par(mar=mar.cj)


plot(conf_alpha$minus_0,ylim=yrange_cj,col = BLUE, type = 'b',main="",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab,cex.axis=cex.axis,las=2)
title(ylab="Confidence",xlab= paste("Within phase groups of",Nphase_trial/Nphase_block,"trials"), line = linelab, cex.lab = cex.lab)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis)
lines(conf_alpha$plus_0, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
lines(conf_alpha$plus_1, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 1)
lines(conf_alpha$minus_1, type = 'b', pch = 16, col = BLUE, lwd = 2, lty = 1)
error.bar(1:length(conf_alpha$minus_0),conf_alpha$minus_0,conf_alpha_sd$minus_0,
          lwd=2, col = BLUE)
error.bar(1:length(conf_alpha$plus_0),conf_alpha$plus_0,conf_alpha_sd$plus_0,
          lwd=2, col = VERMILLION)
error.bar(1:length(conf_alpha$minus_1),conf_alpha$minus_1,conf_alpha_sd$minus_1,
          lwd=2, col = BLUE)
error.bar(1:length(conf_alpha$plus_1),conf_alpha$plus_1,conf_alpha_sd$plus_1,
          lwd=2, col = VERMILLION)

# Plot confidence rolling mean over the course of the experiment ---------------------------------------------------
par(mar=mar.raw.cj)

Nphase <- 6

# Plot ExpB trace 
plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)
Ntrials_phase <- Ntrials/length(unique(Data$phase))

conf_min <- with(subset(cj_ma,sub %in% minus_first),aggregate(cj*6, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% minus_first),aggregate(cj*6, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% plus_first),aggregate(cj*6, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% plus_first),aggregate(cj*6, by=list(trial,cor),se,na.rm=T))
names(conf_plus_se) <- c("trial","cor","cj")

trials_phase1 <- c(seq(0,Ntrials_phase-1), 
                   seq(0,Ntrials_phase-1) + Ntrials_phase*2, 
                   seq(0,Ntrials_phase-1) + Ntrials_phase*4)

xlen <- dim(conf_plus)[1]/2

conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj
conf_plus_cor <- subset(conf_plus,cor==1)$cj
conf_plus_cor_se <- subset(conf_plus_se,cor==1)$cj


plot(conf_min_cor,bty='n',lty = 2,type='l',col="white",ylim=c(4.6,5.5),las=2,     
     main= NULL,cex.lab = cex.lab,cex.axis=cex.axis,xaxt='n',yaxt='n',xlab='',ylab='')
title(xlab = "Trial", ylab = "Confidence", line = linelab, cex.lab = cex.lab)
axis(1,at=seq(0,Ntrials,Ntrials/6),labels=seq(0,Ntrials,Ntrials/6),cex.axis=cex.axis)
axis(2,at=seq(4.6,5.5,.3),labels=seq(4.6,5.5,.3),cex.axis=cex.axis,las=2)

abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')

for (phase in 1:Nphase) {
  if (phase == 1) {
    trials_phase <- seq(Ntrials_phase*(phase-1)+1,Ntrials_phase*phase)
  } else {
    trials_phase <- seq(Ntrials_phase*(phase-1),Ntrials_phase*phase)
  }
  
  if (phase %% 2 == 1) {
    
    lines(trials_phase,conf_min_cor[trials_phase],col=BLUE,lty=2)
    lines(trials_phase,conf_plus_cor[trials_phase],col=VERMILLION)
    
    polygon(c(trials_phase,trials_phase[Ntrials_phase:1]),
            c((conf_min_cor + conf_min_cor_se)[trials_phase],
              (conf_min_cor - conf_min_cor_se)[trials_phase[Ntrials_phase:1]]),
            border=F,col=col_minus_overlay)
    polygon(c(trials_phase,trials_phase[Ntrials_phase:1]),
            c((conf_plus_cor + conf_plus_cor_se)[trials_phase],
              (conf_plus_cor - conf_plus_cor_se)[trials_phase[Ntrials_phase:1]]),
            border=F,col=col_plus_overlay)
  } else {
    
    lines(trials_phase,conf_min_cor[trials_phase],col=VERMILLION,lty=2)
    lines(trials_phase,conf_plus_cor[trials_phase],col=BLUE)
    
    polygon(c(trials_phase,trials_phase[Ntrials_phase:1]),
            c((conf_min_cor + conf_min_cor_se)[trials_phase],
              (conf_min_cor - conf_min_cor_se)[trials_phase[Ntrials_phase:1]]),
            border=F,col=col_plus_overlay)
    polygon(c(trials_phase,trials_phase[Ntrials_phase:1]),
            c((conf_plus_cor + conf_plus_cor_se)[trials_phase],
              (conf_plus_cor - conf_plus_cor_se)[trials_phase[Ntrials_phase:1]]),
            border=F,col=col_minus_overlay)
  }
}  

legend("bottom",border=F,legend=c("High first","Low first"),lty=c(1,2),horiz=T,
       col="black",bty="n",seg.len=1.5,
       cex=cexleg,title = "Counterbalance order")

dev.off()
par(mar=c(5,4,4,2)+.1)



