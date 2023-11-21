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

# Weight traces
# par_trace <- data.frame(trial=(0:(Ntrials-1))+Nskip,sub=rep(subs,each=Ntrials*2),
#                         condition=rep(c("plus","minus"),each=Ntrials),
#                         alpha=0,beta=0,
#                         alpha_alpha_learn=0,beta_alpha_learn=0,
#                         alpha_beta_learn=0,beta_beta_learn=0)


totlen <- length(conditions)*length(difflevels)*length(subs)
# DDM parameters + hyperparameters + error
par_both_learn_fitcj <- data.frame(bound = NA, drift = NA, ter = NA, eta_a = NA,eta_b = NA,
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
  
  ldc_file_both_learn <- paste0('ldc_nn/trim_skip/both_learn/ldcfit_',subs[s],'.Rdata')
  if (file.exists(ldc_file_both_learn)) {
    load(ldc_file_both_learn)
    par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"a0"] <- ldc.results$optim$bestmem[1]
    par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"b0"] <- ldc.results$optim$bestmem[2]
    par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"eta_a"] <- ldc.results$optim$bestmem[4]
    par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"eta_b"] <- ldc.results$optim$bestmem[5]
    par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"cost_ldc"] <- ldc.results$optim$bestval
  } else {
    next
  }
  results <-
    ldc.nn.fit.w(params=c(mean(par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"a0"]),
                          mean(par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"b0"]),1,
                          mean(par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"eta_a"]),
                          mean(par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"eta_b"])),
                 ddm_params = ddm_params,
                 obs=temp_dat,returnFit = F,eta_sep=T,
                 Nupdate_per_trial=Nupdate_per_trial, binning = binning,
                 dt = dt, sigma = sigma)
  anal_sim_both_learn_cjfit[anal_sim_both_learn_cjfit$sim==i&anal_sim_both_learn_cjfit$sub==subs[s] ,'cj'] <- results$pred
  anal_sim_both_learn_cjfit[anal_sim_both_learn_cjfit$sim==i&anal_sim_both_learn_cjfit$sub==subs[s] ,'alpha'] <- results$trace[,1]
  anal_sim_both_learn_cjfit[anal_sim_both_learn_cjfit$sim==i&anal_sim_both_learn_cjfit$sub==subs[s] ,'beta'] <- results$trace[,2]  
  par_both_learn_fitcj[par_both_learn_fitcj$sub==subs[s],"manip"] <- unique(temp_dat$manip)
}

par_both_learn_fitcj <- par_both_learn_fitcj[complete.cases(par_both_learn_fitcj$a0),]
summary(par_both_learn_fitcj)
hist(par_both_learn_fitcj$eta_a,bty='n',main="2 Learning rates model",xlab="Alpha learning rate")
hist(par_both_learn_fitcj$eta_b,bty='n',main="2 Learning rates model",xlab="Beta learning rate")
hist(par_both_learn_fitcj$b0,bty='n',main="2 Learning rates model",xlab="a0")
hist(par_both_learn_fitcj$a0,bty='n',main="2 Learning rates model",xlab="b0")
hist(par_both_learn_fitcj$cost_ldc,bty='n',main="2 Learning rates model",xlab="Cost")


cj_pred <- with(anal_sim_both_learn,aggregate(cj,by=list(trial,sub),mean))
alpha <- with(anal_sim_both_learn,aggregate(alpha,by=list(trial,sub),mean))
beta <- with(anal_sim_both_learn,aggregate(beta,by=list(trial,sub),mean))
names(cj_pred) <- c("trial","sub","cj_pred_both_learn_fitcj")
names(alpha) <- c("trial","sub","alpha_both_learn_fitcj")
names(beta) <- c("trial","sub","beta_both_learn_fitcj")
anal_sim_both_learn_mean <- merge(merge(cj_pred,alpha),beta)
Data <- Data[,!(names(Data) %in% c ("cj_pred_both_learn_fitcj","alpha_both_learn_fitcj","beta_both_learn_fitcj"))]
Data <- merge(Data,anal_sim_both_learn_mean)
Data_fitcj <- Data
Data_fitcj <- Data_fitcj[complete.cases(Data_fitcj$cj_pred_both_learn_fitcj),]
Data_alpha <- subset(Data_fitcj,manip=='alpha')
Data_beta <- subset(Data_fitcj,manip=='beta')

anal_sim_both_learn_cjfit <- anal_sim_both_learn_cjfit[complete.cases(anal_sim_both_learn_cjfit$cj),]

n <- 25 # Rolling mean window size
n_err <- 25

trial_conf_sub <- with(Data,aggregate(cj,by=list(trial,cor,sub),mean))
names(trial_conf_sub) <- c("trial","cor","sub","cj")


pred_conf_sub_both_learn_fitcj <- with(Data_fitcj,aggregate(cj_pred_both_learn_fitcj,by=list(trial,cor,sub),mean))
names(pred_conf_sub_both_learn_fitcj) <- c("trial","cor","sub","cj")


trials <- data.frame(trial=rep((0:(Ntrials-1))+Nskip,each=2),
                     cor=c(0,1),sub=rep(subs,each=Ntrials*2))

cj_pred_ma_both_learn_fitcj <- merge(pred_conf_sub_both_learn_fitcj,trials,all=T)

ma <- function(x,n,names){
  return(rollapply(x[,names], width=n, FUN=function(x) mean(x, na.rm=TRUE),partial=TRUE, align="center"))
}
for (s in subs) {
  cj_pred_ma_both_learn_fitcj[cj_pred_ma_both_learn_fitcj$sub==s&cj_pred_ma_both_learn_fitcj$cor==0,"cj"] <- ma(subset(cj_pred_ma_both_learn_fitcj,sub==s&cor==0),n_err,"cj")
  cj_pred_ma_both_learn_fitcj[cj_pred_ma_both_learn_fitcj$sub==s&cj_pred_ma_both_learn_fitcj$cor==1,"cj"] <- ma(subset(cj_pred_ma_both_learn_fitcj,sub==s&cor==1),n,"cj")
}
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
par_higheta <- subset(par,eta_a>10000)

plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)

# Plot trace alpha experiment
alpha_trace <- with(Data_alpha,aggregate(alpha_both_learn_fitcj,by=list(sub=sub,trial=trial,condition=condition),mean))
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


beta_trace <- with(Data_alpha,aggregate(beta_both_learn_fitcj,by=list(sub=sub,trial=trial,condition=condition),mean))
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

alpha_trace <- with(Data_beta,aggregate(alpha_both_learn_fitcj,by=list(sub=sub,trial=trial,condition=condition),mean))
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


beta_trace <- with(Data_beta,aggregate(beta_both_learn_fitcj,by=list(sub=sub,trial=trial,condition=condition),mean))
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

