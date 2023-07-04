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
sourceCpp("ldc_train.cpp")
# Load dataset --------------------------------------------------------

train_alpha <- read.csv("train_alpha.csv")
train_beta <- read.csv("train_beta.csv")
train <- fastmerge(train_alpha,train_beta)

train$rt2 <- train$rt + train$RTconf
train$cj_pred <- -99
train$cj_pred_eta0 <- -99
train$cj_pred_per_sub <- -99
train$cj_lab <- train$cj
train$cj <- train$cj/6
binning <- F
subs <- unique(train$sub)
conditions <- sort(unique(train$condition))
difflevels <- sort(unique(train$difflevel)) # Important to sort to match with drift order


# Take example subj dataset -----------------------------------------------


s <- 7
obs <- subset(train, sub==s & condition == "plus")
confRTname="RTconf";diffname="difflevel";respname="resp";
totRTname='rt2';targetname='cj';accname='cor';beta_input=.1;error_type='mse'
Nupdate_per_trial <- 1000
binning <- T; nbin <- 6
sigma <- .1
dt <- .001

go_to("fit")
ddm_file1 <- paste0('ddm/ddmfit_',s,'_plus.Rdata')
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
if (sigma==1) {
  ddm_params1 <- ddm_params1*10
}

ldc_file <- paste0('w0/binning/batch_1000/ldcfit_',s,'.Rdata')
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


# Estimate evidence -------------------------------------------------------

drift <- ddm_params1[2:length(ddm_params1)] # Make it more flexible
bound <- ddm_params1[1]
difficulty <- sort(unique(obs$difflevel))

#' Predict single trial amount of evidence accumulated
obs['evidence'] <- bound

if (!("trial" %in% names(obs)) ) {
  obs$trial <- 1:nrow(obs)
}
#' We repeat each trial several times to have the ANN update parameters
#' several times per trial and also to take into account the stochastic nature
#' of estimating single trial evidence accumulation process 
obs_nn <- obs[rep(seq_len(nrow(obs)), each=Nupdate_per_trial), ]

for (trial in seq(1,dim(obs_nn)[1],Nupdate_per_trial)) {
  #' Post decision drift rate sign depends on accuracy 
  if (obs_nn[trial,accname] %in% c(1,'correct','cor')) {
    obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] <- 
      obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] + 
      DDM_fixed_time(v = drift[difficulty==obs_nn[trial,diffname]],
                     time=obs_nn[trial,confRTname],ntrials=Nupdate_per_trial,s=sigma,dt=dt)[,1]
  }else if (obs_nn[trial,accname] %in% c(-1,0,'error','err')) {
    obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] <- 
      obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] + 
      DDM_fixed_time(v = - drift[difficulty==obs_nn[trial,diffname]],
                     time=obs_nn[trial,confRTname],ntrials=Nupdate_per_trial,s=sigma,dt=dt)[,1]
  }
}


Nupdate_1000 <- 1000
obs_nn_1000 <- obs[rep(seq_len(nrow(obs)), each=Nupdate_1000), ]

for (trial in seq(1,dim(obs_nn_1000)[1],Nupdate_1000)) {
  #' Post decision drift rate sign depends on accuracy 
  if (obs_nn_1000[trial,accname] %in% c(1,'correct','cor')) {
    obs_nn_1000[trial:(trial+Nupdate_1000-1),'evidence'] <- 
      obs_nn_1000[trial:(trial+Nupdate_1000-1),'evidence'] + 
      DDM_fixed_time(v = drift[difficulty==obs_nn_1000[trial,diffname]],
                     time=obs_nn_1000[trial,confRTname],ntrials=Nupdate_1000,s=sigma,dt=dt)[,1]
  }else if (obs_nn_1000[trial,accname] %in% c(-1,0,'error','err')) {
    obs_nn_1000[trial:(trial+Nupdate_1000-1),'evidence'] <- 
      obs_nn_1000[trial:(trial+Nupdate_1000-1),'evidence'] + 
      DDM_fixed_time(v = - drift[difficulty==obs_nn_1000[trial,diffname]],
                     time=obs_nn_1000[trial,confRTname],ntrials=Nupdate_1000,s=sigma,dt=dt)[,1]
  }
}
# Train -------------------------------------------------------------------
fac <- 1
#' Input data (evidence and intercept)
x = matrix(c(obs_nn$evidence*fac, #ev
             rep(beta_input,dim(obs_nn)[1])*fac, #bias
             1/sqrt(obs_nn[,totRTname])/fac),  ncol=3) #time

#' Output (confidence)
y = obs_nn[,targetname] 
# y <- y/6

#' x1000 = matrix(c(obs_nn_1000$evidence*fac, #ev
#'                  rep(beta_input,dim(obs_nn_1000)[1])*fac, #bias
#'                  1/sqrt(obs_nn_1000[,totRTname])/fac),  ncol=3) #time
#' 
#' #' Output (confidence)
#' y1000 = obs_nn_1000[,targetname] 
#' 
#' params=c(ldc.results$optim$bestmem[1:3],1)
#' #' Initialize weights
#' w <- params[1:3]
#' w <- c(0,0,1)


# # Best fitting learning rate ----------------------------------------------
# 
# res <- c()
# res_mean <- c()
# res_mode <- c()
# res1000 <- c()
# res_mean1000 <- c()
# res_mode1000 <- c()
# etas <- seq(0,10,.25)
# for (i in etas) {
#   results <- train_model(x,w,y,eta=i,error_type = error_type,trace=F,
#                          binning=T,nbin=nbin,Nupdate_per_trial = Nupdate_per_trial)
#   res <- c(res,results$err)
#   results <- train_model(x1000,w,y1000,eta=i,error_type = error_type,trace=F,
#                          binning=T,nbin=nbin,Nupdate_per_trial = Nupdate_1000)
#   res1000 <- c(res1000,results$err)
#   results <- train_model(x,w,y,eta=i,error_type = error_type,trace=F,cost="per_trial_mean",
#                          binning=T,nbin=nbin,Nupdate_per_trial = Nupdate_per_trial)
#   res_mean <- c(res_mean,results$err)
#   results <- train_model(x1000,w,y1000,eta=i,error_type = error_type,trace=F,cost="per_trial_mean",
#                          binning=T,nbin=nbin,Nupdate_per_trial = Nupdate_1000)
#   res_mean1000 <- c(res_mean1000,results$err)
#   results <- train_model(x,w,y,eta=i,error_type = error_type,trace=F,cost="per_trial_mode",
#                          binning=T,nbin=nbin,Nupdate_per_trial = Nupdate_per_trial)
#   res_mode <- c(res_mode,results$err)
#   results <- train_model(x1000,w,y1000,eta=i,error_type = error_type,trace=F,cost="per_trial_mode",
#                          binning=T,nbin=nbin,Nupdate_per_trial = Nupdate_1000)
#   res_mode1000 <- c(res_mode1000,results$err)
# }
# 
# go_to("plots")
# go_to("check")
# 
# jpeg(filename=paste0("cost_trad_error",s,".jpg"),width=12,height = 8,units="cm",res=300)
# plot(res,type = "p",xaxt='n',xlab='Learning rate',ylab = "Cost", main = "Traditional error")
# axis(1,at=seq(1,length(etas),8),labels = seq(0,10,2))
# points(x=which.min(res),y=min(res),col='red',pch=16)
# points(res1000,pch=4)
# points(x=which.min(res1000),y=min(res1000),col='red',pch=4)
# legend("topright", pch = c(1,4),bty='n',
#        legend=c(paste(Nupdate_per_trial,"updates per trial"),paste(Nupdate_1000,"updates per trial")))
# dev.off()
# 
# jpeg(filename=paste0("cost_agg_mean",s,".jpg"),width=12,height = 8,units="cm",res=300)
# plot(res_mean,type = "p",xaxt='n',xlab='Learning rate',ylab = "Cost", main = "Aggregate trial sims with mean")
# axis(1,at=seq(1,length(etas),8),labels = seq(0,10,2))
# points(x=which.min(res_mean),y=min(res_mean),col='red',pch=16)
# points(res_mean1000,pch=4)
# points(x=which.min(res_mean1000),y=min(res_mean1000),col='red',pch=4)
# legend("topright", pch = c(1,4),bty='n',
#        legend=c(paste(Nupdate_per_trial,"updates per trial"),paste(Nupdate_1000,"updates per trial")))
# dev.off()
# 
# jpeg(filename=paste0("cost_agg_mode",s,".jpg"),width=12,height = 8,units="cm",res=300)
# plot(res_mode,type = "p",xaxt='n',xlab='Learning rate',ylab = "Cost", main = "Aggregate trial sims with mode")
# axis(1,at=seq(1,length(etas),8),labels = seq(0,10,2))
# points(x=which.min(res_mode),y=min(res_mode),col='red',pch=16)
# points(res_mode1000,pch=4)
# points(x=which.min(res_mode1000),y=min(res_mode1000),col='red',pch=4)
# legend("topright", pch = c(1,4),bty='n',
#        legend=c(paste(Nupdate_per_trial,"updates per trial"),paste(Nupdate_1000,"updates per trial")))
# dev.off()
# 
# 
# 
# # Model fits --------------------------------------------------------------
# 
# obs$cj_pred <-
#   ldc.nn.fit.w(params=c(0,0,1,1),ddm_params = ddm_params1,
#                obs=obs,returnFit = F,
#                Nupdate_per_trial=Nupdate_per_trial, binning = binning,
#                dt = dt, sigma = sigma)
# obs$cj_pred2 <- -99
# obs$cj_pred2 <-
#   ldc.nn.fit.w(params=c(0,0,1,1),ddm_params = ddm_params1,
#                obs=obs,returnFit = F,aggreg_pred="mode",
#                Nupdate_per_trial=Nupdate_per_trial, binning = binning,
#                dt = dt, sigma = sigma)
# obs$cj_pred3 <- -99
# obs$cj_pred3 <-
#   ldc.nn.fit.w(params=c(0,0,1,1),ddm_params = ddm_params1,
#                obs=obs,returnFit = F,
#                Nupdate_per_trial=Nupdate_1000, binning = binning,
#                dt = dt, sigma = sigma)
# obs$cj_pred4 <- -99
# obs$cj_pred4 <-
#   ldc.nn.fit.w(params=c(0,0,1,6),ddm_params = ddm_params1,
#                obs=obs,returnFit = F,
#                Nupdate_per_trial=Nupdate_1000, binning = binning,
#                dt = dt, sigma = sigma)
# obs$cj_pred5 <- -99
# obs$cj_pred5 <-
#   ldc.nn.fit.w(params=c(0,0,1,1),ddm_params = ddm_params1,
#                obs=obs,returnFit = F,aggreg_pred="mode",
#                Nupdate_per_trial=Nupdate_1000, binning = binning,
#                dt = dt, sigma = sigma)
# 
# jpeg(filename = paste0("modelfit_mean_mode",s,".jpg"),width=12,height=8,unit='cm',res=300)
# plot(obs$cj_pred,type='l',col='red',ylim=c(0,1),ylab = "Confidence",xlab="Trial",
#      bty='n',main="Model fits") 
# lines(obs$cj,col='black')
# lines(obs$cj_pred2,col='blue')
# legend("bottomleft",lty = c(NA,1,1),pch=c(4,NA,NA),col=c("black","red","blue"),
#        legend = c("Empirical data","Predicted: mean","Predicted: mode"),bty='n',
#        cex=.6)
# dev.off()
# 
# jpeg(filename = paste0("modelfit_learning_rate",s,".jpg"),width=12,height=8,unit='cm',res=300)
# plot(obs$cj_pred4,type='l',col='blue',ylim=c(0,1),ylab = "Confidence",xlab="Trial",
#      bty='n',main="Model fits") 
# points(obs$cj,col='black',pch=4,cex=.5)
# lines(obs$cj_pred3,col='red')
# legend("bottomleft",lty = c(NA,1,1),pch=c(4,NA,NA),col=c("black","red","blue"),
#        legend = c("Empirical data","Predicted: eta = 1","Predicted: eta = 6"),bty='n',
#        cex=.6)
# dev.off()
# 
# 
# jpeg(filename = paste0("modelfit_mean_mode_1000_1",s,".jpg"),width=12,height=8,unit='cm',res=300)
# plot(obs$cj_pred5,type='l',col='blue',ylim=c(0,1),ylab = "Confidence",xlab="Trial",
#      bty='n',main="Model fits") 
# points(obs$cj,col='black',pch=4,cex=.5)
# lines(obs$cj_pred3,col='red')
# legend("bottomleft",lty = c(NA,1,1),pch=c(4,NA,NA),col=c("black","red","blue"),
#        legend = c("Empirical data","Predicted: mean","Predicted: mode"),bty='n',
#        cex=.6)
# dev.off()
# 
# Alpha/Beta traces -------------------------------------------------------

# results1000 <- train_model(x1000,w,y1000,eta=6,error_type = error_type,trace=T,
#                        binning=T,nbin=nbin,Nupdate_per_trial = Nupdate_1000)
results <- train_model(x,ldc.results$optim$bestmem[1:3],y,eta=ldc.results$optim$bestmem[4],error_type = error_type,trace=T,
                       binning=T,nbin=nbin,Nupdate_per_trial = Nupdate_per_trial)
# results1000_1 <- train_model(x1000,w,y1000,eta=1,error_type = error_type,trace=T,
#                            binning=T,nbin=nbin,Nupdate_per_trial = Nupdate_1000)

jpeg(filename = paste0("alpha_trace",s,".jpg"),width=12,height=8,unit='cm',res=300)
plot(results$trace[seq(Nupdate_per_trial,nrow(results$trace),Nupdate_per_trial),1],type='l',
     xlab = "Trial",ylab="Alpha",main = paste("Alpha trace sub",s))
# lines(results$trace[seq(Nupdate_per_trial,nrow(results$trace),Nupdate_per_trial),1],col="orange")
# lines(results1000_1$trace[seq(Nupdate_1000,nrow(results1000_1$trace),Nupdate_1000),1],col="red")
# legend("topright",lty=rep(1,3),col=c("blue","orange","red"),bty='n',cex=.8,
#        legend=c("1000 updates, eta = 6","6000 updates, eta = 1","1000 updates, eta = 1"))
dev.off()

jpeg(filename = paste0("beta_trace",s,".jpg"),width=12,height=8,unit='cm',res=300)
plot(results$trace[seq(Nupdate_per_trial,nrow(results$trace),Nupdate_per_trial),2],type='l',
     xlab = "Trial",ylab="Beta",main = paste("Beta trace sub",s))
# lines(results$trace[seq(Nupdate_per_trial,nrow(results$trace),Nupdate_per_trial),2],col="orange")
# lines(results1000_1$trace[seq(Nupdate_1000,nrow(results1000_1$trace),Nupdate_1000),2],col="red")
# legend("topleft",lty=rep(1,3),col=c("blue","orange","red"),bty='n',cex=.8,
       # legend=c("1000 updates, eta = 6","6000 updates, eta = 1","1000 updates, eta = 1"))
dev.off()

# Don't remember ----------------------------------------------------------
# 
#   for (i in 1:100) {
#     pred <- combine_input_1trial(x[i,],w,binning=T)
#     grad(x[i,], y[i]/6, pred, error_type)
#   }
#   (y[i]/6 - pred) * x[i,1] / sqrt(x[i,3])
#   pred * (1 - pred) * (y[i] - pred) * x[i,1] / sqrt(x[i,3])
#   
#   plot(results$trace[,1])
#   plot(results$trace[,2])
#   
#   a <- 1
#   b <- 1
#   ( 1 /(1 + exp(- x[i,3]* (a*x[i,1] + b/10))))
#   combine_input_1trial(x[i,],c(a,b,1),binning=F)
#   
#   w <- c(a,b,1)
#   z = w[3] * x[i,3] * (w[1] * x[i,1] + w[2] * x[i,2]);
#   z =  1 / (1 + exp(-z));
#   