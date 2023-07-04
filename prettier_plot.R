# Plot traces -------------------------------------------------------------
se <- function(x,na.rm=F) sd(x,na.rm=na.rm)/sqrt(length(x))

width <- 16 # Plot size expressed in cm
height <- 10

go_to("plots")
go_to("w0")
# Plot Exp2A trace --------------------------------------------------------
conf_min <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Alpha_data.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,bty='n',lty = 2,type='l',col='blue',ylim=c(.5,.9),
     main= expression(paste("Exp ",alpha," empirical confidence")),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,lty = 2,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
legend("bottomleft",legend = c(expression(paste("high ",alpha," feedback")), expression(paste("low ",alpha," feedback"))),
       col = c("red","blue"), bty = 'n', lty = c(1,1))
legend("bottomright",legend = c("Correct trials", "Error trials"),
       bty = 'n', lty = c(1,2))
dev.off()
# Plot Exp2A prediction ---------------------------------------------------
conf_min <- with(subset(cj_pred_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma,sub %in% unique(train_alpha$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma,sub %in% unique(train_alpha$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Alpha_pred.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,bty='n',lty = 2,type='l',col='blue',ylim=c(.5,.9),
     main= expression(paste("Exp ",alpha," predicted confidence")),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,lty = 2,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
dev.off()

# Plot Exp2B trace --------------------------------------------------------
conf_min <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Beta_data.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,bty='n',lty = 2,type='l',col='blue',ylim=c(.5,.9),     
     main= expression(paste("Exp ",beta," empirical confidence")),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,lty = 2,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
legend("bottomleft",legend = c(expression(paste("high ",beta," feedback")), expression(paste("low ",beta," feedback"))),
       col = c("red","blue"), bty = 'n', lty = c(1,1))
dev.off()
# Plot Exp2B prediction ---------------------------------------------------
conf_min <- with(subset(cj_pred_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),mean,na.rm=T))
names(conf_min) <- c("trial","cor","cj")
conf_min_se <- with(subset(cj_pred_ma,sub %in% unique(train_beta$sub)),aggregate(minus, by=list(trial,cor),se,na.rm=T))
names(conf_min_se) <- c("trial","cor","cj")
conf_plus <- with(subset(cj_pred_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),mean,na.rm=T))
names(conf_plus) <- c("trial","cor","cj")
conf_plus_se <- with(subset(cj_pred_ma,sub %in% unique(train_beta$sub)),aggregate(plus, by=list(trial,cor),se,na.rm=T))
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

jpeg("Beta_pred.jpg",width=width,height=height,units = 'cm',res=300)
plot(conf_min_err,bty='n',lty = 2,type='l',col='blue',ylim=c(.5,.9),     
     main= expression(paste("Exp ",beta," predicted confidence")),
     xlab = "Trial", ylab = "Confidence")
polygon(c(1:xlen,xlen:1),c(conf_min_err + conf_min_err_se,(conf_min_err - conf_min_err_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_min_cor,col='blue')
polygon(c(1:xlen,xlen:1),c(conf_min_cor + conf_min_cor_se,(conf_min_cor - conf_min_cor_se)[xlen:1]),
        border=F,col=rgb(0,0,255,51,maxColorValue = 255))
lines(conf_plus_err,lty = 2,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_err + conf_plus_err_se,(conf_plus_err - conf_plus_err_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
lines(conf_plus_cor,col='red')
polygon(c(1:xlen,xlen:1),c(conf_plus_cor + conf_plus_cor_se,(conf_plus_cor - conf_plus_cor_se)[xlen:1]),
        border=F,col=rgb(255,0,0,51,maxColorValue = 255))
dev.off()
