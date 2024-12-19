Nphase_block <- 100
Ntrials_phase <- ntrial_per_switch*length(difflevels)
simDat_both$withinphasetrial <- simDat_both$trial %% Ntrials_phase 
simDat_both$phase_block <-  simDat_both$withinphasetrial %/% (Ntrials_phase/Nphase_block)
simDat_both$phase_block <- as.factor(simDat_both$phase_block)


# Aggregating behavior
conf_group <- with(simDat_both,aggregate(cj,by=list(phase_block,manip,condition,cor),mean))
names(conf_group) <- c('phase_block','manip','condition','cor','cj')

trials_phase <- data.frame(phase_block=rep(0:(Nphase_block-1),each=4),
                           cor=c(0,1),condition=rep(c("minus","plus"),each=2),manip=rep(fb_manip,each=Nphase_block*4))

conf_group <- merge(conf_group,trials_phase,all=T)
table(complete.cases(conf_group$cj))
count_complete <- function(dat) {
  return(sum(complete.cases(dat)))
}

conf_alpha_count <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor,count_complete)
conf_alpha <- cast(subset(conf_group,manip=='alpha'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,fun.aggregate = mean,na.rm=T)
conf_beta_count <- cast(subset(conf_group,manip=='beta'),phase_block~condition+cor,count_complete)


plot(conf_alpha$minus_0,ylim=c(.5,.9),col = BLUE, type = 'b',main="",
     lty = 2, pch = 17, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",cex.main=cex.title,
     xlab = paste("Consecutive groups of",Ntrials_phase/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
lines(conf_alpha$plus_0, type = 'b', pch = 17, col = VERMILLION, lwd = 2, lty = 2)
lines(conf_alpha$plus_1, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 1)
lines(conf_alpha$minus_1, type = 'b', pch = 16, col = BLUE, lwd = 2, lty = 1)
legend("top",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)

plot(conf_beta$minus_0,ylim=c(.4,.9),col = BLUE, type = 'b',main="",
     lty = 2, pch = 17, lwd = 2, bty = 'n', xaxt = 'n', ylab = "Confidence",cex.main=cex.title,
     xlab = paste("Consecutive groups of",Ntrials_phase/Nphase_block,"trials"),cex.lab=cex.lab*.66/.83,cex.axis=cex.axis*.66/.83)
axis(1, at = 1:Nphase_block, labels = 1:Nphase_block,cex.axis=cex.axis*.66/.83)
lines(conf_beta$plus_0, type = 'b', pch = 17, col = VERMILLION, lwd = 2, lty = 2)
lines(conf_beta$plus_1, type = 'b', pch = 16, col = VERMILLION, lwd = 2, lty = 1)
lines(conf_beta$minus_1, type = 'b', pch = 16, col = BLUE, lwd = 2, lty = 1)
legend("top",legend = c('High','Low'),lty = c(1,1),col = c(VERMILLION,BLUE),
       pch = c(16,16),horiz = T, bty = 'n',cex = cex.legend*.66/.83)

