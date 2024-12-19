
# Plot model prediction + correlation with empirical per subject ----------
# Correlation between empirical and model-predicted confidence per subject
grouped_alpha <- split(Data_alpha[,c('cj','cj_pred_both_learn')],list(Data_alpha$sub))
correlations_alpha <- sapply(grouped_alpha,function(x) cor(x$cj,x$cj_pred_both_learn))
correlations_alpha <- correlations_alpha[!is.na(correlations_alpha)]
grouped_beta <- split(Data_beta[,c('cj','cj_pred_both_learn')],list(Data_beta$sub))
correlations_beta <- sapply(grouped_beta,function(x) cor(x$cj,x$cj_pred_both_learn))
correlations_beta <- correlations_beta[!is.na(correlations_beta)]
cor_range <- range(c(correlations_alpha,correlations_beta))


jpeg(filename = paste0("traces_poster_alpha.jpg"),units = 'cm',width = 18,height = 19,res=300)
# layout(matrix(c(1,1,3,3,2,2,4,4,5,6,7,8),ncol=3))
layout(matrix(c(1,3,2,4),ncol=1),heights = c(.05,.4,.05,.4))

par(mar=c(0,0,0,0))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Model Fits"))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Goodness of fit"))
par(mar=c(4,5,0,0)+.1)


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


hist(correlations_alpha,main="",xlab="Empirical - Model correlation",
     ylab="Nb subjects",xlim=cor_range,cex.lab = cex.lab,cex.axis=cex.axis)
# Add a thick black vertical line at the mean
abline(v=mean(correlations_alpha),lwd=2,col="black")

dev.off()

jpeg(filename = paste0("traces_poster_beta.jpg"),units = 'cm',width = 18,height = 19,res=300)
# layout(matrix(c(1,1,3,3,2,2,4,4,5,6,7,8),ncol=3))
layout(matrix(c(1,3,2,4),ncol=1),heights = c(.05,.4,.05,.4))

par(mar=c(0,0,0,0))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Model Fits"))
plot.new()
title(cex.main=cex.title,line=title_line,main= expression("Goodness of fit"))
par(mar=c(4,5,0,0)+.1)

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

hist(correlations_beta,main="",xlab="Empirical - Model correlation",
     ylab="Nb subjects",xlim=cor_range,cex.lab = cex.lab,cex.axis=cex.axis)
# Add a thick black vertical line at the mean
abline(v=mean(correlations_beta),lwd=2,col="black")
dev.off()


# Plot aggregated traces  -------------------------------------------------
m <- 'both'
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


xlen <- nrow(conf_alpha)
go_to("plots")
tiff(paste0('trace_aggreg_poster_alpha.tiff'), width = 18, height = 26, units = 'cm', res = 300)
layout(matrix(c(1,2,3),ncol=1),heights = c(1,.5,.5))

par(mar=c(0,5,4,2)+.1)

plot(conf_alpha$minus_0,ylim=c(.7,.9),col = BLUE, type = 'b',main=expression(paste(alpha,"-Manipulated Feedback")),
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",cex.main=cex.title,
     # xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"),
     xlab = '',
     cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, 
     # labels = 1:Nphase_block,
     labels = F,
     cex.axis=cex.axis*.66/.83)
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

par(mar=c(0,5,4,2)+.1)

plot(rt_alpha$minus,ylim=c(.8,1.2),col = BLUE, type = 'b',main="",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, 
     # labels = 1:Nphase_block,
     labels=F,
     cex.axis=cex.axis*.66/.83
     )
title(ylab = "Reaction time (s)", 
      # xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"), 
      line = 2.5,cex.lab=cex.lab*.66/.83)
lines(rt_alpha$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
error.bar(1:xlen,rt_alpha$minus,rt_alpha_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:xlen,rt_alpha$plus,rt_alpha_sd$plus,
          lwd=2, col = VERMILLION)

par(mar=c(4,5,0,2)+.1)

plot(cor_alpha$minus,ylim=c(.5,1),col = BLUE, type = 'b',main="",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
title(ylab = "Accuracy", xlab = paste("Within phase groups of",Nphase_trial/Nphase_block,"trials"), 
      line = 2.5,cex.lab=cex.lab*.66/.83)
lines(cor_alpha$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
error.bar(1:xlen,cor_alpha$minus,cor_alpha_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:xlen,cor_alpha$plus,cor_alpha_sd$plus,
          lwd=2, col = VERMILLION)

dev.off()

tiff(paste0('trace_aggreg_poster_beta.tiff'), width = 18, height = 26, units = 'cm', res = 300)
layout(matrix(c(1,2,3),ncol=1),heights = c(1,.5,.5))
par(mar=c(0,5,4,2)+.1)

plot(conf_beta$minus_0,ylim=c(.7,.9),col = BLUE, type = 'b',main=expression(paste(beta,"-Manipulated Feedback")),
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = F,cex.axis=cex.axis*.66/.83)
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
       pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)

par(mar=c(0,5,4,2)+.1)

plot(rt_beta$minus,ylim=c(.8,1.2),col = BLUE, type = 'b',main="",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = F,cex.axis=cex.axis*.66/.83)
title(ylab = "Reaction time (s)", xlab = paste("Consecutive groups of",Nphase_trial/Nphase_block,"trials"), 
      line = 2.5,cex.lab=cex.lab*.66/.83)
lines(rt_beta$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
error.bar(1:xlen,rt_beta$minus,rt_beta_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:xlen,rt_beta$plus,rt_beta_sd$plus,
          lwd=2, col = VERMILLION)

par(mar=c(4,5,0,2)+.1)

plot(cor_beta$minus,ylim=c(.5,1),col = BLUE, type = 'b',main="",
     lty = 2, pch = 16, lwd = 2, bty = 'n', xaxt = 'n', ylab = "",cex.main=cex.title,
     xlab = "",cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
title(ylab = "Accuracy", xlab = paste("Within phase groups of",Nphase_trial/Nphase_block,"trials"), 
      line = 2.5,cex.lab=cex.lab*.66/.83)
lines(cor_beta$plus, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 2)
error.bar(1:xlen,cor_beta$minus,cor_beta_sd$minus,
          lwd=2, col = BLUE)
error.bar(1:xlen,cor_beta$plus,cor_beta_sd$plus,
          lwd=2, col = VERMILLION)

dev.off()

