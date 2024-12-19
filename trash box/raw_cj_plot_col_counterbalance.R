# Plot confidence rolling mean over the course of the experiment ---------------------------------------------------
par(mar=mar.raw.cj)

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

xlen <- dim(conf_plus)[1]/2
conf_min_err <- subset(conf_min,cor==0)$cj
conf_min_err_se <- subset(conf_min_se,cor==0)$cj
conf_min_cor <- subset(conf_min,cor==1)$cj
conf_min_cor_se <- subset(conf_min_se,cor==1)$cj
conf_plus_err <- subset(conf_plus,cor==0)$cj
conf_plus_err_se <- subset(conf_plus_se,cor==0)$cj
conf_plus_cor <- subset(conf_plus,cor==1)$cj
conf_plus_cor_se <- subset(conf_plus_se,cor==1)$cj

plot(conf_min_err,bty='n',lty = 2,type='l',col=BLUISH_GREEN,ylim=c(3.6,5.4),las=2,     
     main= NULL,cex.lab = cex.lab,cex.axis=cex.axis,xaxt='n',yaxt='n',xlab='',ylab='')
title(xlab = "Trial", ylab = "Confidence", line = linelab, cex.lab = cex.lab)
axis(1,at=seq(0,Ntrials,Ntrials/6),labels=seq(0,Ntrials,Ntrials/6),cex.axis=cex.axis)
axis(2,at=seq(3.6,5.4,.3),labels=seq(3.6,5.4,.3),cex.axis=cex.axis,las=2)
abline(v=seq(Ntrials_phase,Ntrials-1,Ntrials_phase),lty=2,col='lightgrey')
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,158,115,51,maxColorValue = 255))
lines(conf_min_cor,col=BLUISH_GREEN)
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,158,115,51,maxColorValue = 255))
lines(conf_plus_err,lty = 2,col=REDDISH_PURPLE)
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(204,121,167,51,maxColorValue = 255))
lines(conf_plus_cor,col=REDDISH_PURPLE)
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(204,121,167,51,maxColorValue = 255))
legend(cex=cexleg,"bottomleft",legend = c(paste("High feedback first"), paste("Low feedback first")),
       col = c(REDDISH_PURPLE,BLUISH_GREEN), bty = 'n', lty = c(1,1))
legend("bottomright",border=F,legend=c("Correct","Incorrect"),lty=c(1,2),horiz=F,
       col="black",bty="n",seg.len=1.5,
       cex=cexleg,title = "Trial accuracy")

dev.off()
par(mar=c(5,4,4,2)+.1)

