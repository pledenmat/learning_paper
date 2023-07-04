# Investigate alpha/beta traces -------------------------------------------
par$alpha <- NA
par$beta <- NA
for (s in subs) {
  for (cond in conditions) {
    print(paste("Retrieving alpha/beta for sub",s,cond,"condition"))
    temp_dat <- subset(train,sub==s&condition==cond)
    temp_par <- subset(par,sub==s&condition==cond)
    
    drift <- temp_par$drift 
    difficulty <- temp_par$difflevel
    
    temp_dat$evidence <- mean(temp_par$bound)
    
    temp_dat_nn <- temp_dat[rep(seq_len(nrow(temp_dat)), each=Nupdate_per_trial), ]
    
    for (trial in seq(1,dim(temp_dat_nn)[1],Nupdate_per_trial)) {
      #' Post decision drift rate sign depends on accuracy 
      if (temp_dat_nn[trial,"cor"] %in% c(1,'correct','cor')) {
        temp_dat_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] <- 
          temp_dat_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] + 
          DDM_fixed_time(v = drift[difficulty==temp_dat_nn[trial,"difflevel"]],
                         time=temp_dat_nn[trial,"RTconf"],ntrials=Nupdate_per_trial,s=sigma,dt=dt)[,1]
      }else if (temp_dat_nn[trial,"cor"] %in% c(-1,0,'error','err')) {
        temp_dat_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] <- 
          temp_dat_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] + 
          DDM_fixed_time(v = - drift[difficulty==temp_dat_nn[trial,"difflevel"]],
                         time=temp_dat_nn[trial,"RTconf"],ntrials=Nupdate_per_trial,s=sigma,dt=dt)[,1]
      }
    }
    
    #' Input data (evidence and intercept)
    x = matrix(c(temp_dat_nn$evidence, #ev
                 rep(beta_input,dim(temp_dat_nn)[1]), #bias
                 1/sqrt(temp_dat_nn$rt2)),  ncol=3) #time
    
    #' Output (confidence)
    y = temp_dat_nn$cj 
    
    #' Initialize weights
    w <- w0
    
    results <- train_model(x,w,y,eta=mean(temp_par$eta),error_type = "mse",trace=F)
    
    par[par$sub==s&par$condition==cond,"alpha"] <- results$w[1]
    par[par$sub==s&par$condition==cond,"beta"] <- results$w[2]
  }
}

# Look at estimated parameters -----------------------------------------
go_to("plots")

# Learning rate
jpeg(paste0("eta_",Nupdate_per_trial,".jpg"),width=20,height=15,units='c <- ',res=300)
hist(par$eta, xlab = "Learning rate", 
     main = paste("Estimated learning rate,",Nupdate_per_trial,"updates per trial"))
dev.off()

breaks_alpha <- seq(-2,14,2)
breaks_beta <- seq(-10,30,5)
# Alpha
jpeg("Alpha_hist_exp2a.jpg",width=20,height=15,units='cm',res=300)
hist(subset(par,condition=="minus"&manip=="alpha")$alpha,col=rgb(1,0,0,.5), xlab = "Alpha", 
     main = "Estimated alpha - Alpha experiment",breaks=breaks_alpha)
hist(subset(par,condition=="plus"&manip=="alpha")$alpha,add=T, col = rgb(0,0,1,.5),breaks=breaks_alpha)
legend("topright",bty = "n",fill = c(rgb(0,0,1,.5),rgb(1,0,0,.5)),legend = c("Minus","Plus"))
dev.off()

jpeg("Alpha_hist_exp2b.jpg",width=20,height=15,units='cm',res=300)
hist(subset(par,condition=="minus"&manip=="beta")$alpha,col=rgb(1,0,0,.5), xlab = "Alpha", 
     main = "Estimated alpha - Beta experiment",breaks=breaks_alpha)
hist(subset(par,condition=="plus"&manip=="beta")$alpha,add=T, col = rgb(0,0,1,.5),breaks=breaks_alpha)
legend("topright",bty = "n",fill = c(rgb(0,0,1,.5),rgb(1,0,0,.5)),legend = c("Minus","Plus"))
dev.off()

# Beta
jpeg("Beta_hist_exp2a.jpg",width=20,height=15,units='cm',res=300)
hist(subset(par,condition=="minus"&manip=="alpha")$beta,col=rgb(1,0,0,.5), xlab = "Beta", 
     main = "Estimated beta - Alpha experiment",breaks=breaks_beta)
hist(subset(par,condition=="plus"&manip=="alpha")$beta,add=T, col = rgb(0,0,1,.5),breaks=breaks_beta)
legend("topright",bty = "n",fill = c(rgb(0,0,1,.5),rgb(1,0,0,.5)),legend = c("Minus","Plus"))
dev.off()

jpeg("Beta_hist_exp2b.jpg",width=20,height=15,units='cm',res=300)
hist(subset(par,condition=="minus"&manip=="beta")$beta,col=rgb(1,0,0,.5), xlab = "Beta", 
     main = "Estimated beta - Beta experiment",breaks=breaks_beta)
hist(subset(par,condition=="plus"&manip=="beta")$beta,add=T, col = rgb(0,0,1,.5),breaks=breaks_beta)
legend("topright",bty = "n",fill = c(rgb(0,0,1,.5),rgb(1,0,0,.5)),legend = c("Minus","Plus"))
dev.off()

m <- lmer(data = par, alpha ~ condition*manip + (1|sub))
anova(m)
m <- lmer(data = par, beta ~ condition*manip + (1|sub))
anova(m)

m <- lmer(data = subset(par,manip=="alpha"), alpha ~ condition + (1|sub))
anova(m)
m <- lmer(data = subset(par,manip=="alpha"), beta ~ condition + (1|sub))
anova(m)

m <- lmer(data = subset(par,manip=="beta"), alpha ~ condition + (1|sub))
anova(m)
emm <- emmeans(m)
m <- lmer(data = subset(par,manip=="beta"), beta ~ condition + (1|sub))
anova(m)

