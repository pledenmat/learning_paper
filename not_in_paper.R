# Not in paper -------------------------------------------------------------------
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

conf_min <- with(subset(cj_pred_ma_skip,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_skip,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_skip,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_skip,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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

conf_min <- with(subset(cj_pred_ma_skip,sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma_skip,sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma_skip,sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma_skip,sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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
# Plot traces LR fixed to 5 ---------------------------------------------------

width <- 16 # Plot size expressed in cm
height <- 10

se <- function(x,na.rm=F) sd(x,na.rm=na.rm)/sqrt(length(x))
beta_range <- c(5,15)
alpha_range <- c(5,15)
title_line <- -2
cex.title <- 3
cex.lab <- 3
cex.axis <- 2
cex.legend <- 2
go_to("plots")
go_to("alternating_fb")
go_to("trim")
for (m in models) {
  jpeg(filename = paste0("traces_",m,"_learn_5.jpg"),units = 'cm',width = 42,height = 30,res=300)
  # layout(matrix(c(1,1,3,3,2,2,4,4,5,6,7,8),ncol=3))
  layout(matrix(c(1,2,9,9,5,6,11,11,1,3,10,10,5,7,12,12,1,4,13,14,5,8,15,16),ncol=3),heights = c(.05,.05,.2,.2,.05,.05,.2,.2))
  par(mar=c(0,0,0,0))
  plot.new()
  title(cex.main=cex.title,line=title_line,main = expression(paste(alpha,"-Manipulated Feedback, LR = 5")))
  plot.new()
  title(cex.main=cex.title,line=title_line,main= expression("Empirical Data"))
  plot.new()
  title(cex.main=cex.title,line=title_line,main= expression("Model Fits"))
  plot.new()
  title(cex.main=cex.title,line=title_line,main= expression("Weight traces"))
  plot.new()
  title(cex.main=cex.title,line=title_line,main = expression(paste(beta,"-Manipulated Feedback, LR = 5")))
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
  
  conf_min <- with(subset(cj_pred_fixed_lr_5,model==m&sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
  names(conf_min) <- c("trial","cor","cj")
  conf_min_se <- with(subset(cj_pred_fixed_lr_5,model==m&sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
  names(conf_min_se) <- c("trial","cor","cj")
  conf_plus <- with(subset(cj_pred_fixed_lr_5,model==m&sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
  names(conf_plus) <- c("trial","cor","cj")
  conf_plus_se <- with(subset(cj_pred_fixed_lr_5,model==m&sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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
  
  conf_min <- with(subset(cj_pred_fixed_lr_5,model==m&sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
  names(conf_min) <- c("trial","cor","cj")
  conf_min_se <- with(subset(cj_pred_fixed_lr_5,model==m&sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
  names(conf_min_se) <- c("trial","cor","cj")
  conf_plus <- with(subset(cj_pred_fixed_lr_5,model==m&sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
  names(conf_plus) <- c("trial","cor","cj")
  conf_plus_se <- with(subset(cj_pred_fixed_lr_5,model==m&sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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
  alpha_trace <- Data_alpha[,c(paste0("alpha_",m,"_learn_fixed_lr_5"),"sub","trial","condition")]
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
  
  
  beta_trace <- Data_alpha[,c(paste0("beta_",m,"_learn_fixed_lr_5"),"sub","trial","condition")]
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
  
  alpha_trace <- Data_beta[,c(paste0("alpha_",m,"_learn_fixed_lr_5"),"sub","trial","condition")]
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
  
  
  beta_trace <- Data_beta[,c(paste0("beta_",m,"_learn_fixed_lr_5"),"sub","trial","condition")]
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

# Plot traces LR fixed to 1 ---------------------------------------------------

width <- 16 # Plot size expressed in cm
height <- 10

se <- function(x,na.rm=F) sd(x,na.rm=na.rm)/sqrt(length(x))
beta_range <- c(5,15)
alpha_range <- c(5,15)
title_line <- -2
cex.title <- 3
cex.lab <- 3
cex.axis <- 2
cex.legend <- 2
go_to("plots")
go_to("alternating_fb")
go_to("trim")
for (m in models) {
  jpeg(filename = paste0("traces_",m,"_learn_1.jpg"),units = 'cm',width = 42,height = 30,res=300)
  # layout(matrix(c(1,1,3,3,2,2,4,4,5,6,7,8),ncol=3))
  layout(matrix(c(1,2,9,9,5,6,11,11,1,3,10,10,5,7,12,12,1,4,13,14,5,8,15,16),ncol=3),heights = c(.05,.05,.2,.2,.05,.05,.2,.2))
  par(mar=c(0,0,0,0))
  plot.new()
  title(cex.main=cex.title,line=title_line,main = expression(paste(alpha,"-Manipulated Feedback, LR = 1")))
  plot.new()
  title(cex.main=cex.title,line=title_line,main= expression("Empirical Data"))
  plot.new()
  title(cex.main=cex.title,line=title_line,main= expression("Model Fits"))
  plot.new()
  title(cex.main=cex.title,line=title_line,main= expression("Weight traces"))
  plot.new()
  title(cex.main=cex.title,line=title_line,main = expression(paste(beta,"-Manipulated Feedback, LR = 1")))
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
  
  conf_min <- with(subset(cj_pred_fixed_lr_1,model==m&sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
  names(conf_min) <- c("trial","cor","cj")
  conf_min_se <- with(subset(cj_pred_fixed_lr_1,model==m&sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
  names(conf_min_se) <- c("trial","cor","cj")
  conf_plus <- with(subset(cj_pred_fixed_lr_1,model==m&sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
  names(conf_plus) <- c("trial","cor","cj")
  conf_plus_se <- with(subset(cj_pred_fixed_lr_1,model==m&sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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
  
  conf_min <- with(subset(cj_pred_fixed_lr_1,model==m&sub %in% minus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
  names(conf_min) <- c("trial","cor","cj")
  conf_min_se <- with(subset(cj_pred_fixed_lr_1,model==m&sub %in% minus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
  names(conf_min_se) <- c("trial","cor","cj")
  conf_plus <- with(subset(cj_pred_fixed_lr_1,model==m&sub %in% plus_first),aggregate(cj, by=list(trial,cor),mean,na.rm=T))
  names(conf_plus) <- c("trial","cor","cj")
  conf_plus_se <- with(subset(cj_pred_fixed_lr_1,model==m&sub %in% plus_first),aggregate(cj, by=list(trial,cor),se,na.rm=T))
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
  alpha_trace <- Data_alpha[,c(paste0("alpha_",m,"_learn_fixed_lr_1"),"sub","trial","condition")]
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
  
  
  beta_trace <- Data_alpha[,c(paste0("beta_",m,"_learn_fixed_lr_1"),"sub","trial","condition")]
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
  
  alpha_trace <- Data_beta[,c(paste0("alpha_",m,"_learn_fixed_lr_1"),"sub","trial","condition")]
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
  
  
  beta_trace <- Data_beta[,c(paste0("beta_",m,"_learn_fixed_lr_1"),"sub","trial","condition")]
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


# One big experiment ------------------------------------------------------

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


