# Plot traces -------------------------------------------------------------
se <- function(x,na.rm=F) sd(x,na.rm=na.rm)/sqrt(length(x))

title_line <- -2
cex.title <- 3
cex.lab <- 3
cex.axis <- 2
cex.legend <- 2
go_to("plots")
go_to("alternating_fb")
go_to("trim_skip")

# jpeg(filename = "traces_cor.jpg",units = 'cm',width = 39,height = 30,res=300)
tiff(filename = "traces_cor_2col.tiff",units = 'cm',width = 36,height = 30,res=300)
# OG
# layout(matrix(c(1,2,9,9,5,6,11,11,1,3,10,10,5,7,12,12,1,4,13,14,5,8,15,16),ncol=3),heights = c(.05,.05,.2,.2,.05,.05,.2,.2))
# Cor
# layout(matrix(c(1,2,9,5,6,11,1,3,10,5,7,12,1,4,13,5,8,14),ncol=3),heights = c(.05,.05,.4,.05,.05,.4))
# Cor 2 col
layout(matrix(c(1,2,9,3,10,4,13,5,6,11,7,12,8,14),ncol=2),heights = c(.05,.05,.4,.05,.4,.05,.4))
# Only trace
# layout(matrix(c(1,2,7,4,5,9,1,3,8,4,6,10),ncol=2),heights = c(.05,.05,.4,.05,.05,.4))
par(mar=c(0,0,0,0))
plot.new()
title(cex.main=cex.title,line=title_line,main = expression(paste(alpha,"-Manipulated Feedback")))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Empirical Data"))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Model Fits"))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Goodness of fit"))
plot.new()
title(cex.main=cex.title,line=title_line,main = expression(paste(beta,"-Manipulated Feedback")))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Empirical Data"))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Model Fits"))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Goodness of fit"))
par(mar=c(4,5,0,0)+.1)
# Plot ExpA trace --------------------------------------------------------

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

# isna <- function(data) {
#   return(sum(is.na(data)))
# }
# test <- with(cj_pred_ma,aggregate(cj, by=list(sub,cor),isna))
# names(test) <- c("sub","cor","count")
# test <- cast(test, sub ~ cor)
# summary(test)
# hist(subset(test,cor==1)$count)
# hist(subset(test,cor==0)$count)
xlen <- dim(conf_plus)[1]/2
conf_min_err <- subset(conf_min,cor==0)$cj
conf_min_err_se <- subset(conf_min_se,cor==0)$cj
conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj
conf_plus_err <- subset(conf_plus,cor==0)$cj
conf_plus_err_se <- subset(conf_plus_se,cor==0)$cj
conf_plus_cor <- subset(conf_plus,cor==1)$cj
conf_plus_cor_se <- subset(conf_plus_se,cor==1)$cj

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.6,.9),
     main= NULL,cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=round(seq(Ntrials_phase,Ntrials-1,Ntrials_phase)),lty=2,col='black',lwd=2)
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

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.6,.9),
     main= NULL,cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=round(seq(Ntrials_phase,Ntrials-1,Ntrials_phase)),lty=2,col='black',lwd=2)
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

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.6,.9),     
     main= NULL,cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=round(seq(Ntrials_phase,Ntrials-1,Ntrials_phase)),lty=2,col='black',lwd=2)
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
legend(cex=cex.legend,"bottomleft",legend = c(expression(paste("high ",beta," feedback first")), 
                                              expression(paste("low ",beta," feedback first"))),
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

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUE,ylim=c(.6,.9),     
     main= NULL,cex.lab = cex.lab,cex.axis=cex.axis,
     xlab = "Trial", ylab = "Confidence")
abline(v=round(seq(Ntrials_phase,Ntrials-1,Ntrials_phase)),lty=2,col='black',lwd=2)
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

# Plot correlation --------------------------------------------------------
Data_alpha$cj_plot <- jitter(Data_alpha$cj)
## Use densCols() output to get density at each point
x <- densCols(Data_alpha$cj_plot,Data_alpha$cj_pred, colramp=colorRampPalette(c("black", "white")))
Data_alpha$dens <- col2rgb(x)[1,] + 1L

## Map densities to colors
cols <-  grey(seq(.9,0,length.out=256))
Data_alpha$col <- cols[Data_alpha$dens]

## Plot it, reordering rows so that densest points are plotted on top
plot(cj_pred~cj_plot, data=Data_alpha[order(Data_alpha$dens),], pch=20, col=col,
     xaxt='n',bty='n',xlab="Empirical confidence",ylab = "Model confidence",
     cex.axis=cex.axis, cex.lab=cex.lab)
axis(1,at=(1:6)/6,labels = 1:6,cex.axis=cex.axis)
abline(lm(Data_alpha$cj_pred~Data_alpha$cj)$coef,lwd=2)
mtext(side = 3, adj = 1, cex=2,
      paste0('r=',round(cor(Data_alpha$cj,Data_alpha$cj_pred,method='spearman'),3)))

Data_beta$cj_plot <- jitter(Data_beta$cj)
## Use densCols() output to get density at each point
x <- densCols(Data_beta$cj_plot,Data_beta$cj_pred, colramp=colorRampPalette(c("black", "white")))
Data_beta$dens <- col2rgb(x)[1,] + 1L

## Map densities to colors
cols <-  grey(seq(.9,0,length.out=256))
Data_beta$col <- cols[Data_beta$dens]

## Plot it, reordering rows so that densest points are plotted on top
plot(cj_pred~cj_plot, data=Data_beta[order(Data_beta$dens),], pch=20, col=col,
     xaxt='n',bty='n',xlab="Empirical confidence",ylab = "Model confidence",
     cex.axis=cex.axis, cex.lab=cex.lab)
axis(1,at=(1:6)/6,labels = 1:6,cex.axis=cex.axis)
abline(lm(Data_beta$cj_pred~Data_beta$cj)$coef,lwd=2)
mtext(side = 3, adj = 1, cex=2,
      paste0('r=',round(cor(Data_beta$cj,Data_beta$cj_pred,method='spearman'),3)))
dev.off()


# Plot only parameters ----------------------------------------------------


jpeg(filename = "traces_par.jpg",units = 'cm',width = 36,height = 15,res=300)
layout(matrix(c(1,3,4,2,5,6),ncol=2),heights=c(1,4,4))
par(mar=c(0,0,0,0))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Weight traces"))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Weight traces"))
par(mar=c(4,5,0,0)+.1)

# Plot parameter traces ---------------------------------------------------
plus_first <- unique(subset(Data_alpha,phase==0&condition=='plus')$sub)
minus_first <- unique(subset(Data_alpha,phase==0&condition=='minus')$sub)

# Plot trace alpha experiment
alpha_trace <- with(subset(par_trace,sub %in% Data_alpha$sub),aggregate(alpha,by=list(sub=sub,trial=trial,condition=condition),mean))
alpha_trace_minus <- cast(subset(alpha_trace,sub %in% minus_first),sub ~trial, value = "x", fun.aggregate = mean)
alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first ),sub~trial, value = "x", fun.aggregate = mean)
count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
count_plus <- count_plus[2:length(count_plus)]
count_minus <- count_minus[2:length(count_minus)]

plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
     ylim = c(-10,30),bty='n')
abline(v=round(seq(Ntrials_phase,Ntrials-1,Ntrials_phase)),lty=2,col='black',lwd=2)
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
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(0,30),bty='n')
abline(v=round(seq(Ntrials_phase,Ntrials-1,Ntrials_phase)),lty=2,col='black',lwd=2)
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
alpha_trace_plus <- cast(subset(alpha_trace,sub %in% plus_first),sub~trial, value = "x", fun.aggregate = mean)
count_plus <- sapply(alpha_trace_plus, function(y) sum(length(which(!is.na(y)))))
count_minus <- sapply(alpha_trace_minus, function(y) sum(length(which(!is.na(y)))))
count_plus <- count_plus[2:length(count_plus)]
count_minus <- count_minus[2:length(count_minus)]

plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(alpha_trace_minus,na.rm=T),type='l',col=BLUE,xlab="",ylab="Alpha",
     ylim=c(-10,25),bty='n')
abline(v=round(seq(Ntrials_phase,Ntrials-1,Ntrials_phase)),lty=2,col='black',lwd=2)
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
beta_trace_plus <- cast(subset(beta_trace,sub %in% plus_first),sub~trial, value = "x", fun.aggregate = mean)
plot(cex.lab = cex.lab,cex.axis=cex.axis,colMeans(beta_trace_minus,na.rm=T),type='l',col=BLUE,xlab="Trial",ylab="Beta",
     ylim=c(0,30),bty='n')
abline(v=round(seq(Ntrials_phase,Ntrials-1,Ntrials_phase)),lty=2,col='black',lwd=2)
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





