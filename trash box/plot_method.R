alpha_plus <- rep(c(36,NA),each=126,length.out=756)
alpha_minus <- rep(c(NA,9),each=126,length.out=756)
beta_plus <- rep(c(0,NA),each=126,length.out=756)
beta_minus <- rep(c(NA,0),each=126,length.out=756)
lwdline <- 3

tiff(filename = 'method.tiff',width = 21,height = 13, units = 'cm', res = 300)
par(mfrow=c(1,2),mar=c(5,4,1,1))
plot(alpha_plus,ylim=c(9,36),xlim=c(0,800),type='l',bty='n', lwd=lwdline,
     ylab='',xlab='',col=VERMILLION,cex.axis=cexax)
mtext(side=1,line=2.5,text='Trial',cex=cexlab)
mtext(side=2,line=2.5,text='Alpha',cex=cexlab)
lines(alpha_minus,lwd=lwdline,col=BLUE)
plot(beta_plus,xlim=c(0,800),type='l',bty='n', lwd=lwdline,
     ylab='',xlab='',col=VERMILLION,cex.axis=cexax,cex.lab=cexlab)
mtext(side=1,line=2.5,text='Trial',cex=cexlab)
mtext(side=2,line=2.5,text='Beta',cex=cexlab)
lines(beta_minus,lwd=lwdline,col=BLUE)
dev.off()
par(mfrow=c(1,1),mar=c(5,4,4,2)+.1)
