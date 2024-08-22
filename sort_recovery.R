
# Recovery among subjects who had learning model as best one --------------

for (m in models) {
  temp_par <- subset(par,model==m&sub %in% sub_learn_best)
  temp_par <- with(temp_par,aggregate(cbind(eta_a,eta_b,a0,b0),list(sub=sub),mean))
  temp_fit_par <- subset(fit_par,gen_model==m& fit_model==m &sub %in% sub_learn_best)
  # reorder both according to subject
  temp_fit_par <- temp_fit_par[order(temp_fit_par$sub),]
  temp_par <- temp_par[order(temp_par$sub),]
  par(mfrow=c(2,2))
  plot(temp_par$eta_a,temp_fit_par$eta_a,xlab="Generative",ylab="Estimated",
       main=paste("eta_a, cor =",
                  round(cor(method = "spearman",temp_par$eta_a,temp_fit_par$eta_a),3),"p =",
                  round(cor.test(method = "spearman",temp_par$eta_a,temp_fit_par$eta_a)$p.value,3)))
  plot(temp_par$eta_b,temp_fit_par$eta_b,xlab="Generative",ylab="Estimated",
       main=paste("eta_b, cor =",
                  round(cor(method = "spearman",temp_par$eta_b,temp_fit_par$eta_b),3),"p =",
                  round(cor.test(method = "spearman",temp_par$eta_b,temp_fit_par$eta_b)$p.value,3)))
  plot(temp_par$a0,temp_fit_par$a0,xlab="Generative",ylab="Estimated",
       main=paste("Alpha 0, cor =",
                  round(cor(method = "spearman",temp_par$a0,temp_fit_par$a0),3),"p =",
                  round(cor.test(method = "spearman",temp_par$a0,temp_fit_par$a0)$p.value,3)))
  plot(temp_par$b0,temp_fit_par$b0,xlab="Generative",ylab="Estimated",
       main=paste("Beta 0, cor =",
                  round(cor(method = "spearman",temp_par$b0,temp_fit_par$b0),3),"p =",
                  round(cor.test(method = "spearman",temp_par$b0,temp_fit_par$b0)$p.value,3)))
  
}
par(mfrow=c(1,1))


# Model comparison only dynamic group -----------------------------------
fit_par$bic <- bic_custom(fit_par$cost_ldc,fit_par$Npar,fit_par$Ndata_point)
mean_bic <- with(subset(fit_par,sub %in% sub_learn_best),aggregate(bic,by=list(fit=fit_model,gen=gen_model),mean))
mean_bic <- cast(mean_bic,fit~gen)

bic_sub <- with(subset(fit_par,sub %in% sub_learn_best),aggregate(bic,by=list(fit=fit_model,gen=gen_model,sub=sub),mean))
bic_sub <- cast(bic_sub,gen+sub~fit)
bic_sub$win_model <- sort(models)[apply(bic_sub[,3:6],1,which.min)]
with(bic_sub,aggregate(win_model,list(gen=gen),table))
table(subset(bic_sub,gen=="no")$win_model)
table(subset(bic_sub,gen=="alpha")$win_model)
table(subset(bic_sub,gen=="beta")$win_model)
table(subset(bic_sub,gen=="both")$win_model)


