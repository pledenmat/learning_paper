# Retrieve fits -----------------------------------------------------------
# Some fixed arguments
beta_input <- .1
Nupdate_per_trial <- 1
model <- "allpar"
dt <- .001; sigma <- .1
binning <- F
Nsim <- 1

conditions <- sort(unique(Data$condition))
difflevels <- sort(unique(Data$difflevel)) # Important to sort to match with drift order

totlen <- length(conditions)*length(difflevels)*length(subs)
# DDM parameters + hyperparameters + error
par_both_learn_fitcj <- data.frame(bound = NA, drift = NA, ter = NA, eta_a = NA,eta_b = NA,
                             cost_ddm = NA, cost_ldc = NA, a0 = NA, b0 = NA,
                             sub = rep(subs, each = length(conditions)*length(difflevels)),
                             condition = rep(conditions,each = length(difflevels),
                                             length.out=totlen),
                             difflevel = rep(difflevels, length.out = totlen),
                             manip=NA)

par_both_learn_wrong_resp <- data.frame(bound = NA, drift = NA, ter = NA, eta_a = NA,eta_b = NA,
                                        cost_ddm = NA, cost_ldc = NA, a0 = NA, b0 = NA,
                                        sub = rep(subs, each = length(conditions)*length(difflevels)),
                                        condition = rep(conditions,each = length(difflevels),
                                                        length.out=totlen),
                                        difflevel = rep(difflevels, length.out = totlen),
                                        manip=NA)
  
go_to("fit")
go_to("alternating_fb")
if (!file.exists("anal_sim_both_learn_cjfit.Rdata")) {
  anal_sim_both_learn_cjfit <- Data[rep(seq_len(nrow(Data)), each=Nsim), c('trial','withinphasetrial','sub','condition','cor')]
  anal_sim_both_learn_cjfit$alpha <- NA
  anal_sim_both_learn_cjfit$beta <- NA
  anal_sim_both_learn_cjfit$cj <- NA
  anal_sim_both_learn_cjfit$sim <- 1:Nsim
}
if (!file.exists("anal_sim_both_learn_wrong_resp.Rdata")) {
  anal_sim_both_learn_wrong_resp <- Data[rep(seq_len(nrow(Data)), each=Nsim), c('trial','withinphasetrial','sub','condition','cor')]
  anal_sim_both_learn_wrong_resp$alpha <- NA
  anal_sim_both_learn_wrong_resp$beta <- NA
  anal_sim_both_learn_wrong_resp$cj <- NA
  anal_sim_both_learn_wrong_resp$sim <- 1:Nsim
}
for (s in 1:length(subs)) {
  print(paste("Retrieving participant",s))
  
  temp_dat <- subset(Data,sub==subs[s])
  
  ddm_file <- paste0('ddm/ddmfit_',subs[s],'.Rdata')
  if(file.exists(ddm_file)){
    load(ddm_file)
  }else{
    next
  }
  ddm_params <- ddm.results$optim$bestmem[c(1,5:length(ddm.results$optim$bestmem))]
  
  ldc_file_both_learn <- paste0('ldc_nn/trim/both_learn/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file_both_learn)) {
    load(ldc_file_both_learn)
    par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
    par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
    par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
    par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
    par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  } else {
    if (s == 56) {
      par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"a0"] <- 2.038468
      par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"b0"] <- 9.438675
      par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"eta_a"] <- 0.002523
      par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"eta_b"] <- 0.001841
      par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"cost_ldc"] <- 0.011340
    } else if (s == 6) {
      par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"a0"] <- 1.457321
      par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"b0"] <- 13.830877
      par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"eta_a"] <- 2.486172
      par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"eta_b"] <- 0.008269
      par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"cost_ldc"] <- 0.007257
    }
  }
  
  ldc_file_both_learn <- paste0('ldc_nn/trim/both_learn/wrong_resp_coding/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file_both_learn)) {
    load(ldc_file_both_learn)
    par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
    par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
    par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
    par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
    par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  } else {
    if (s == 55) {
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"a0"] <- 1.155483
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"b0"] <- 7.369901
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"eta_a"] <- 3.084824
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"eta_b"] <- 1.024432
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"cost_ldc"] <- 0.009181
    } 
    else if (s == 6) {
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"a0"] <- 3.595170
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"b0"] <- 13.683037
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"eta_a"] <- 2.193099
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"eta_b"] <- 0.040854
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"cost_ldc"] <- 0.007276
    }
    else if (s == 7) {
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"a0"] <- 2.955426
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"b0"] <- 8.801150
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"eta_a"] <- 0.008193
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"eta_b"] <- 0.291868
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"cost_ldc"] <- 0.005520
    } 
    else if (s == 10) {
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"a0"] <- 8.100477
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"b0"] <- 13.206168
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"eta_a"] <- 0.003759
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"eta_b"] <- 0.007621
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"cost_ldc"] <- 0.005640
    }
    else if (s == 23) {
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"a0"] <- 3.989513
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"b0"] <- 9.398379
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"eta_a"] <- 1.042615
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"eta_b"] <- 0.1069040
      par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"cost_ldc"] <- 0.009072
    }
  }
  for (sim in 1:Nsim) {
    if (!file.exists("anal_sim_both_learn_cjfit.Rdata")) {
      results <-
        ldc.nn.fit.w(params=c(mean(par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"a0"]),
                              mean(par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"b0"]),1,
                              mean(par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"eta_a"]),
                              mean(par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"eta_b"])),
                     ddm_params = ddm_params,
                     obs=temp_dat,returnFit = F,eta_sep=T,
                     Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                     dt = dt, sigma = sigma)
      anal_sim_both_learn_cjfit[anal_sim_both_learn_cjfit$sim==sim&anal_sim_both_learn_cjfit$sub==subs[s] ,'cj'] <- results$pred
      anal_sim_both_learn_cjfit[anal_sim_both_learn_cjfit$sim==sim&anal_sim_both_learn_cjfit$sub==subs[s] ,'alpha'] <- results$trace[,1]
      anal_sim_both_learn_cjfit[anal_sim_both_learn_cjfit$sim==sim&anal_sim_both_learn_cjfit$sub==subs[s] ,'beta'] <- results$trace[,2]  
    }
    if (!file.exists("anal_sim_both_learn_wrong_resp.Rdata")) {
      results <-
        ldc.nn.fit.w(params=c(mean(par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"a0"]),
                              mean(par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"b0"]),1,
                              mean(par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"eta_a"]),
                              mean(par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"eta_b"])),
                     ddm_params = ddm_params,
                     obs=temp_dat,returnFit = F,eta_sep=T,
                     Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                     dt = dt, sigma = sigma)
      anal_sim_both_learn_wrong_resp[anal_sim_both_learn_wrong_resp$sim==sim&anal_sim_both_learn_wrong_resp$sub==subs[s] ,'cj'] <- results$pred
      anal_sim_both_learn_wrong_resp[anal_sim_both_learn_wrong_resp$sim==sim&anal_sim_both_learn_wrong_resp$sub==subs[s] ,'alpha'] <- results$trace[,1]
      anal_sim_both_learn_wrong_resp[anal_sim_both_learn_wrong_resp$sim==sim&anal_sim_both_learn_wrong_resp$sub==subs[s] ,'beta'] <- results$trace[,2]  
    }
    
  }
  
  par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  par_both_learn_wrong_resp[par_both_learn_wrong_resp$sub==subs[s],"manip"] <- unique(temp_dat$manip)
  
}
if (!file.exists("anal_sim_both_learn_wrong_resp.Rdata")) {
  save(anal_sim_both_learn_wrong_resp,file="anal_sim_both_learn_wrong_resp.Rdata")
}
if (!file.exists("anal_sim_both_learn_cjfit.Rdata")) {
  save(anal_sim_both_learn_cjfit,file="anal_sim_both_learn_cjfit.Rdata")
}
load("anal_sim_both_learn_cjfit.Rdata")
load("anal_sim_both_learn_wrong_resp.Rdata")

cj_pred <- with(anal_sim_both_learn_cjfit,aggregate(cj,by=list(trial,sub),mean))
alpha <- with(anal_sim_both_learn_cjfit,aggregate(alpha,by=list(trial,sub),mean))
beta <- with(anal_sim_both_learn_cjfit,aggregate(beta,by=list(trial,sub),mean))
names(cj_pred) <- c("trial","sub","cj_pred_both_learn_fitcj")
names(alpha) <- c("trial","sub","alpha_both_learn_fitcj")
names(beta) <- c("trial","sub","beta_both_learn_fitcj")
anal_sim_both_learn_cjfit_mean <- merge(merge(cj_pred,alpha),beta)
Data <- Data[,!(names(Data) %in% c ("cj_pred_both_learn_fitcj","alpha_both_learn_fitcj","beta_both_learn_fitcj"))]
Data <- merge(Data,anal_sim_both_learn_cjfit_mean)

cj_pred <- with(anal_sim_both_learn_wrong_resp,aggregate(cj,by=list(trial,sub),mean))
alpha <- with(anal_sim_both_learn_wrong_resp,aggregate(alpha,by=list(trial,sub),mean))
beta <- with(anal_sim_both_learn_wrong_resp,aggregate(beta,by=list(trial,sub),mean))
names(cj_pred) <- c("trial","sub","cj_pred_both_learn_wrong_resp")
names(alpha) <- c("trial","sub","alpha_both_learn_wrong_resp")
names(beta) <- c("trial","sub","beta_both_learn_wrong_resp")
anal_sim_both_learn_wrong_resp_mean <- merge(merge(cj_pred,alpha),beta)
Data <- Data[,!(names(Data) %in% c ("cj_pred_both_learn_wrong_resp","alpha_both_learn_wrong_resp","beta_both_learn_wrong_resp"))]
Data <- merge(Data,anal_sim_both_learn_wrong_resp_mean)

Data_alpha <- subset(Data,manip=='alpha')
Data_beta <- subset(Data,manip=='beta')

n <- 25 # Rolling mean window size
n_err <- 25

trial_conf_sub <- with(Data,aggregate(cj,by=list(trial,cor,sub),mean))
names(trial_conf_sub) <- c("trial","cor","sub","cj")


pred_conf_sub_both_learn_fitcj <- with(Data,aggregate(cj_pred_both_learn_fitcj,by=list(trial,cor,sub),mean))
names(pred_conf_sub_both_learn_fitcj) <- c("trial","cor","sub","cj")

pred_conf_sub_both_learn_wrong_resp <- with(Data,aggregate(cj_pred_both_learn_wrong_resp,by=list(trial,cor,sub),mean))
names(pred_conf_sub_both_learn_wrong_resp) <- c("trial","cor","sub","cj")

trials <- data.frame(trial=rep((0:(Ntrials-1))+Nskip,each=2),
                     cor=c(0,1),sub=rep(subs,each=Ntrials*2))

cj_ma <- merge(trial_conf_sub,trials,all=T)
cj_pred_ma_both_learn_fitcj <- merge(pred_conf_sub_both_learn_fitcj,trials,all=T)
cj_pred_ma_both_learn_wrong_resp <- merge(pred_conf_sub_both_learn_wrong_resp,trials,all=T)

ma <- function(x,n,names){
  return(rollapply(x[,names], width=n, FUN=function(x) mean(x, na.rm=TRUE),partial=TRUE, align="center"))
}
for (s in subs) {
  cj_ma[cj_ma$sub==s&cj_ma$cor==0,"cj"] <- ma(subset(cj_ma,sub==s&cor==0),n_err,"cj")
  cj_ma[cj_ma$sub==s&cj_ma$cor==1,"cj"] <- ma(subset(cj_ma,sub==s&cor==1),n,"cj")
  cj_pred_ma_both_learn_fitcj[cj_pred_ma_both_learn_fitcj$sub==s&cj_pred_ma_both_learn_fitcj$cor==0,"cj"] <- ma(subset(cj_pred_ma_both_learn_fitcj,sub==s&cor==0),n_err,"cj")
  cj_pred_ma_both_learn_fitcj[cj_pred_ma_both_learn_fitcj$sub==s&cj_pred_ma_both_learn_fitcj$cor==1,"cj"] <- ma(subset(cj_pred_ma_both_learn_fitcj,sub==s&cor==1),n,"cj")
  cj_pred_ma_both_learn_wrong_resp[cj_pred_ma_both_learn_wrong_resp$sub==s&cj_pred_ma_both_learn_wrong_resp$cor==0,"cj"] <- ma(subset(cj_pred_ma_both_learn_wrong_resp,sub==s&cor==0),n_err,"cj")
  cj_pred_ma_both_learn_wrong_resp[cj_pred_ma_both_learn_wrong_resp$sub==s&cj_pred_ma_both_learn_wrong_resp$cor==1,"cj"] <- ma(subset(cj_pred_ma_both_learn_wrong_resp,sub==s&cor==1),n,"cj")
  
}
# Compare both fits ------------------------------------------------------
bic_custom <- function(Residuals,k,n){
  return(log(n)*k+n*log(Residuals/n))
}
par_both_learn_fitcj$model <- "fitcj"
par_both_learn_wrong_resp$model <- "wrong_resp"
par <- rbind(par_both_learn_fitcj,par_both_learn_wrong_resp)
par$Npar <- 4
par$Ndata_point <- nrow(Data)/Nsub

par$bic <- bic_custom(par$cost_ldc,par$Npar,par$Ndata_point)
mean_bic <- with(par,aggregate(bic,by=list(model=model,manip=manip),mean))
mean_bic$delta <- -99
mean_bic[mean_bic$manip=="alpha","delta"] <- 
  mean_bic[mean_bic$manip=="alpha",]$x - 
  min(mean_bic[mean_bic$manip=="alpha",]$x)
mean_bic[mean_bic$manip=="beta","delta"] <- 
  mean_bic[mean_bic$manip=="beta",]$x -
  min(mean_bic[mean_bic$manip=="beta",]$x)

plot(par_both_learn_fitcj$a0 ~ par_both_learn_wrong_resp$a0, bty='n', xlab = "Fit 1 (wrong response coding)", ylab = "Fit 2",
     main = paste("a0, r =",round(cor(par_both_learn_fitcj$a0, par_both_learn_wrong_resp$a0,method="spearman"),3)))
abline(coef=c(0,1))

plot(par_both_learn_fitcj$b0 ~ par_both_learn_wrong_resp$b0, bty='n', xlab = "Fit 1 (wrong response coding)", ylab = "Fit 2",
     main = paste("b0, r =",round(cor(par_both_learn_fitcj$b0, par_both_learn_wrong_resp$b0,method="spearman"),3)))
abline(coef=c(0,1))

plot(par_both_learn_fitcj$eta_a ~ par_both_learn_wrong_resp$eta_a, bty='n', xlab = "Fit 1 (wrong response coding)", ylab = "Fit 2",
     main = paste("eta_a, r =",round(cor(par_both_learn_fitcj$eta_a, par_both_learn_wrong_resp$eta_a,method="spearman"),3)))
abline(coef=c(0,1))

plot(par_both_learn_fitcj$eta_b ~ par_both_learn_wrong_resp$eta_b, bty='n', xlab = "Fit 1 (wrong response coding)", ylab = "Fit 2",
     main = paste("eta_b, r =",round(cor(par_both_learn_fitcj$eta_b, par_both_learn_wrong_resp$eta_b,method="spearman"),3)))
abline(coef=c(0,1))
plot(par_both_learn_fitcj$cost_ldc ~ par_both_learn_wrong_resp$cost_ldc, bty='n', xlab = "Fit 1 (wrong response coding)", ylab = "Fit 2",
     main = paste("cost_ldc, r =",round(cor(par_both_learn_fitcj$cost_ldc, par_both_learn_wrong_resp$cost_ldc,method="spearman"),3)))
abline(coef=c(0,1))

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

jpeg(filename = "traces_both_learn_fitcj_fitcj.jpg",units = 'cm',width = 42,height = 30,res=300)
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

conf_min <- with(subset(cj_pred_ma_both_learn_fitcj,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_both_learn_fitcj,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_both_learn_fitcj,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_both_learn_fitcj,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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

conf_min <- with(subset(cj_pred_ma_both_learn_fitcj,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_both_learn_fitcj,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_both_learn_fitcj,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_both_learn_fitcj,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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
alpha_trace <- with(Data_alpha,aggregate(alpha_both_learn_fitcj,by=list(sub=sub,trial=trial,condition=condition),mean))
alpha_trace_minus <- cast(subset(alpha_trace,sub %in% minus_first),sub ~trial, value = "x", fun.aggregate = mean)
alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first ),sub~trial, value = "x", fun.aggregate = mean)
count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
count_plus <- count_plus[2:length(count_plus)]
count_minus <- count_minus[2:length(count_minus)]

plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
     ylim = c(0,20),bty='n')
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


beta_trace <- with(Data_alpha,aggregate(beta_both_learn_fitcj,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(5,15),bty='n')
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

alpha_trace <- with(Data_beta,aggregate(alpha_both_learn_fitcj,by=list(sub=sub,trial=trial,condition=condition),mean))
alpha_trace_minus <- cast(subset(alpha_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first),sub~trial, value = "x", fun.aggregate = mean)
count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
count_plus <- count_plus[2:length(count_plus)]
count_minus <- count_minus[2:length(count_minus)]

plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
     ylim=c(0,20),bty='n')
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


beta_trace <- with(Data_beta,aggregate(beta_both_learn_fitcj,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(5,15),bty='n')
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


# Plot traces both learn - Wrong response coding --------------------------------------------------
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

jpeg(filename = "traces_both_learn_wrong_resp.jpg",units = 'cm',width = 42,height = 30,res=300)
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

conf_min <- with(subset(cj_pred_ma_both_learn_wrong_resp,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_both_learn_wrong_resp,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_both_learn_wrong_resp,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_both_learn_wrong_resp,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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

conf_min <- with(subset(cj_pred_ma_both_learn_wrong_resp,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_both_learn_wrong_resp,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_both_learn_wrong_resp,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_both_learn_wrong_resp,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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
alpha_trace <- with(Data_alpha,aggregate(alpha_both_learn_wrong_resp,by=list(sub=sub,trial=trial,condition=condition),mean))
alpha_trace_minus <- cast(subset(alpha_trace,sub %in% minus_first),sub ~trial, value = "x", fun.aggregate = mean)
alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first ),sub~trial, value = "x", fun.aggregate = mean)
count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
count_plus <- count_plus[2:length(count_plus)]
count_minus <- count_minus[2:length(count_minus)]

plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
     ylim = c(0,20),bty='n')
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


beta_trace <- with(Data_alpha,aggregate(beta_both_learn_wrong_resp,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(5,15),bty='n')
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

alpha_trace <- with(Data_beta,aggregate(alpha_both_learn_wrong_resp,by=list(sub=sub,trial=trial,condition=condition),mean))
alpha_trace_minus <- cast(subset(alpha_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first),sub~trial, value = "x", fun.aggregate = mean)
count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
count_plus <- count_plus[2:length(count_plus)]
count_minus <- count_minus[2:length(count_minus)]

plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
     ylim=c(0,20),bty='n')
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


beta_trace <- with(Data_beta,aggregate(beta_both_learn_wrong_resp,by=list(sub=sub,trial=trial,condition=condition),mean))
beta_trace_minus <- cast(subset(beta_trace,sub %in% minus_first),sub~trial, value = "x", fun.aggregate = mean)
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(5,15),bty='n')
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


hist(par_both_learn_fitcj$eta_a)
hist(par_both_learn_fitcj$eta_b)
hist(par_both_learn_fitcj$a0)
hist(par_both_learn_fitcj$b0)

# Diagnostic --------------------------------------------------------------
beta_input <- .1
Nupdate_per_trial <- 1
Nsim_err <- 1000
# bound, ter, z, vratio, drifts
dt <- .001; sigma <- .1
error_type1 <- "cross-entropy"
error_type2 <- "mse"
target <- "fb"
fitname <- 'cj'

s <- Nsub
fitted_par <- colMeans(subset(par_both_learn_fitcj,sub==subs[s])[,c("a0","b0","eta_a","eta_b")])
fitted_par <- c(fitted_par[1:2],1,fitted_par[3:4])
og_cost <- mean(subset(par_both_learn_fitcj,sub==subs[s])$cost_ldc)
ddm_file <- paste0('ddm/ddmfit_',subs[s],'.Rdata')
if(file.exists(ddm_file)){ load(ddm_file) }
ddm_params <- ddm.results$optim$bestmem[c(1,5:length(ddm.results$optim$bestmem))]
obs <- subset(Data,sub==subs[s])
obs$ev <- as.numeric(str_sub(obs$ev,2,-2))
# First check that the cost is the same when computed by hand
err <- ldc.nn.fit.w(params=fitted_par,obs,ddm_params,dt=.001,sigma=0.1,Nsim_error=1000,
                    Nupdate_per_trial=1,returnFit=T,estimate_evidence = T,
                    confRTname="RTconf",diffname="difflevel",respname="resp",
                    totRTname='rt2',targetname='fb',accname='cor',beta_input=.1,
                    error_type1='cross-entropy',error_type2='mse',binning=F,nbin=6,
                    shuffle=F,cost="separated",aggreg_pred="mean",Nskip_error=0,
                    eta_sep=T, fitname='cj')
err/og_cost
# Now try some other parameterization
params <- c(fitted_par[1:3],0,-10)
params <- fitted_par
err_explo <- ldc.nn.fit.w(params=params,obs,ddm_params,dt=.001,sigma=0.1,Nsim_error=1000,
                    Nupdate_per_trial=1,returnFit=T,estimate_evidence = T,
                    confRTname="RTconf",diffname="difflevel",respname="resp",
                    totRTname='rt2',targetname='fb',accname='cor',beta_input=.1,
                    error_type1='cross-entropy',error_type2='mse',binning=F,nbin=6,
                    shuffle=F,cost="separated",aggreg_pred="mean",Nskip_error=0,
                    eta_sep=T, fitname='cj')
err_explo/err

colmap <- viridis::viridis(10000)
summary(obs_nn$diff)
j <- 0
y_pred = combine_input(x_err[((j*Nsim_error +1):((j+1)*Nsim_error)),], w, binning = binning, nbin = nbin)
y_pred_mean = mean(y_pred)
err_it <- error(y_err[j+1], y_pred_mean, error_type2)
with(Data,aggregate(cj_pred_both_learn_fitcj,list(cor,condition),mean))
with(Data,aggregate(cj,list(cor,condition),mean))

ggplot(obs_nn, aes(x = ev_adj, y = evidence, colour = diff)) +
  geom_point() + scale_color_viridis(limits=c(-1,1))


plot(obs_nn$ev_adj~obs_nn$evidence, )

estimate_evidence = T;
confRTname="RTconf";diffname="difflevel";respname="resp";
totRTname='rt2';targetname='fb';accname='cor';beta_input=.1;
error_type1='cross-entropy';error_type2='mse';binning=F;nbin=6;
shuffle=F;cost="separated";aggreg_pred="mean";Nskip_error=0;
eta_sep=T; fitname='cj';Nsim_error=1000
returnFit=T

with(obs_nn,aggregate(evidence,list(cor,difflevel,response),mean))
with(obs_nn,aggregate(ev,list(cor,difflevel,response),mean))
test <- with(obs_err,aggregate(evidence,list(cor,difflevel,response),mean))
m <- lm(data = obs_err, evidence ~ response*cor*difflevel)
Anova(m,type = 3)

params[1:2] <- c(0,0)
j <- 0
y_pred = combine_input(x_err[(j*Nsim_error+1):((j+1)*Nsim_error),], w, binning = binning, nbin = nbin)
y_pred_mean = mean(y_pred)
err = error(y_err[j+1], y_pred_mean, error_type2);
library(DEoptim)
optimal_params <- DEoptim(ldc.nn.fit.w,ddm_params=ddm_params, 
                          obs = obs,
                          lower = c(0,-100,1,0,0), 
                          upper = c(50,100,1,10000,10000),
                          Nupdate_per_trial=Nupdate_per_trial,
                          dt = dt, sigma = sigma,binning=binning,eta_sep=T,
                          Nsim_err=Nsim_err,error_type1=error_type1,fitname="fb",
                          error_type2=error_type2,targetname=target,cost="separated",
                          control=c(itermax=1000,steptol=70,reltol=.00001,NP=40))
