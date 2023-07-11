library(Rcpp)
sourceCpp("DDM_fixed_time.cpp")

ldc.nn.fit.all <- function(params,obs,ntrials=10,dt=.001,sigma=.1,Nupdate_per_trial=1000,w0=c(10,-1,1),
                       confRTname="RTconf",diffname="difflevel",respname="resp",
                       totRTname='rt2',targetname='cj',accname='cor',beta_input=.1,error_type='mse'){
  #' Step 1 : have some DDM parameters + a/b
  #' Step 2 : Generate simulated data from DDM
  #' Step 3 : Evaluate DDM fit
  #' Step 4 : Use DDM bound and drift rate to infer evidence accumulated at each trial
  #' Step 5 : Gradiant descent 
  #' Step 6 : Global cost is DDM cost + NN cost
  #' Loop and optimize DDM parameters + NN hyperparameters (learning rate + batch size)

  drift <- params[7:length(params)] # Make it more flexible
  params <- params[1:6]
  names(params) <- c('a','ter','z','sigma','vratio','eta')
  difficulty <- sort(unique(obs[,diffname]))
  
  #' First generate DDM predictions
  for (d in 1:length(drift)) {
    if (d == 1) {
      predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(
        v = drift[d], a = params['a'], ter = params['ter'], z = params['z'],
        ntrials = ntrials*dim(obs)[1]/length(drift), s = sigma, dt = dt, postdriftmod = params['vratio'], 
        t2distribution = rep(obs[,confRTname],times=ntrials)
      ))
      names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
      predictions['drift'] <- drift[d]
      predictions[diffname] <- difficulty[d]
    }
    else{
      temp_pred <- data.frame(DDM_with_confidence_slow_fullconfRT(
        v = drift[d], a = params['a'], ter = params['ter'], z = params['z'],
        ntrials = ntrials*dim(obs)[1]/length(drift), s = sigma, dt = dt, postdriftmod = params['vratio'], 
        t2distribution = rep(obs[,confRTname],times=ntrials)
      ))
      names(temp_pred) <- c('rt','resp','cor','evidence2','rt2','cj')
      temp_pred['drift'] <- drift[d]
      temp_pred[diffname] <- difficulty[d]
      predictions <- fastmerge(predictions,temp_pred)
      
    }
  }
  
  #' Evaluate goodness-of-fit with RT quantiles for correct and errors
  c_pred <- predictions[predictions$cor == 1,]
  e_pred <- predictions[predictions$cor == 0,]
  
  c_obs <- obs[obs[,accname] %in% c(1,'cor','correct'),]
  e_obs <- obs[obs[,accname] %in% c(-1,0,'err','error'),]
  
  obs_props <- NULL; pred_props <- NULL
  for (d in 1:length(drift)) {
    #' Now, get the quantile RTs on the "obs data" for correct and error distributions 
    #' separately (for quantiles .1, .3, .5, .7, .9)
    c_quantiles <- quantile(c_obs[c_obs[,diffname] == difficulty[d],]$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    e_quantiles <- quantile(e_obs[e_obs[,diffname] == difficulty[d],]$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    if (any(is.na(e_quantiles))) {
      e_quantiles <- rep(0,5)
    }
    if (any(is.na(c_quantiles))) {
      c_quantiles <- rep(0,5)
    }      
    #' to combine correct and incorrect we scale the expected interquantile 
    #' probability by the proportion of correct and incorrect respectively
    prop_obs_c <- dim(c_obs[c_obs[,diffname] == difficulty[d],])[1] / dim(obs)[1]
    prop_obs_e <- dim(e_obs[e_obs[,diffname] == difficulty[d],])[1] / dim(obs)[1]
    
    c_obs_proportion = prop_obs_c * c(.1, .2, .2, .2, .2, .1)
    e_obs_proportion = prop_obs_e * c(.1, .2, .2, .2, .2, .1)
    obs_props <- c(obs_props,c_obs_proportion,e_obs_proportion)
    
    c_pred_rt <- sort(c_pred[c_pred[,diffname] == difficulty[d],]$rt)
    e_pred_rt <- sort(e_pred[e_pred[,diffname] == difficulty[d],]$rt)
    # now, get the proportion of responses that fall between the obs quantiles when applied to the pred data (scale by N?)
    c_pred_proportion <- c(
      sum(c_pred_rt <= c_quantiles[1]),
      sum(c_pred_rt <= c_quantiles[2]) - sum(c_pred_rt <= c_quantiles[1]),
      sum(c_pred_rt <= c_quantiles[3]) - sum(c_pred_rt <= c_quantiles[2]),
      sum(c_pred_rt <= c_quantiles[4]) - sum(c_pred_rt <= c_quantiles[3]),
      sum(c_pred_rt <= c_quantiles[5]) - sum(c_pred_rt <= c_quantiles[4]),
      sum(c_pred_rt > c_quantiles[5])
    ) / dim(predictions)[1]
    
    e_pred_proportion <- c(
      sum(e_pred_rt <= e_quantiles[1]),
      sum(e_pred_rt <= e_quantiles[2]) - sum(e_pred_rt <= e_quantiles[1]),
      sum(e_pred_rt <= e_quantiles[3]) - sum(e_pred_rt <= e_quantiles[2]),
      sum(e_pred_rt <= e_quantiles[4]) - sum(e_pred_rt <= e_quantiles[3]),
      sum(e_pred_rt <= e_quantiles[5]) - sum(e_pred_rt <= e_quantiles[4]),
      sum(e_pred_rt > e_quantiles[5])
    ) / dim(predictions)[1]
    pred_props <- c(pred_props,c_pred_proportion,e_pred_proportion)
  }
  ddm_cost <- sum((obs_props-pred_props)^2)
  
  #' Predict single trial amount of evidence accumulated
  obs['evidence'] <- params['a']
  #' We repeat each trial several times to have the ANN update parameters
  #' several times per trial and also to take into account the stochastic nature
  #' of estimating single trial evidence accumulation process 
  obs_nn <- obs[rep(seq_len(nrow(obs)), each=Nupdate_per_trial), ]

  for (trial in seq(1,dim(obs_nn)[1],Nupdate_per_trial)) {
    #' Post decision drift rate sign depends on accuracy 
    if (obs_nn[trial,accname] %in% c(1,'correct','cor')) {
      obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] <- 
        obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] + 
        DDM_fixed_time(v = drift[difficulty==obs_nn[trial,diffname]],
                       time=obs_nn[trial,confRTname],ntrials=Nupdate_per_trial,s=sigma,dt=dt)[,1]
    }else if (obs_nn[trial,accname] %in% c(-1,0,'error','err')) {
      obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] <- obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] + 
        DDM_fixed_time(v = - drift[difficulty==obs_nn[trial,diffname]],
                       time=obs_nn[trial,confRTname],ntrials=Nupdate_per_trial,s=sigma,dt=dt)[,1]
    }
  }
  
  #' Input data (evidence and intercept)
  x = matrix(c(obs_nn$evidence, #ev
               rep(beta_input,dim(obs_nn)[1]), #bias
               1/sqrt(obs_nn[,totRTname])),  ncol=3) #time
  
  #' Output (confidence)
  y = obs_nn[,targetname] 
  
  #' Initialize weights
  w <- w0
  
  results <- train_model(x,w,y,eta=params['eta'],error_type = error_type,
                         Nupdate_per_trial = 1,trace=F)
  
  cost <- ddm_cost*10 + results$err
  return(cost)
}

ddm.fit <- function(params,obs,ntrials=10,dt=.001,sigma=.1,
                    confRTname="RTconf",diffname="difflevel",accname='cor'){
  #' Step 1 : have some DDM parameters + a/b
  #' Step 2 : Generate simulated data from DDM
  #' Step 3 : Evaluate DDM fit
  #' Step 4 : Use DDM bound and drift rate to infer evidence accumulated at each trial
  #' Step 5 : Gradiant descent 
  #' Step 6 : Global cost is DDM cost + NN cost
  #' Loop and optimize DDM parameters + NN hyperparameters (learning rate + batch size)
  
  drift <- params[5:length(params)] # Make it more flexible
  params <- params[1:4]
  names(params) <- c('a','ter','z','vratio')
  difficulty <- sort(unique(obs[,diffname]))
  
  #' First generate DDM predictions
  for (d in 1:length(drift)) {
    if (d == 1) {
      predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(
        v = drift[d], a = params['a'], ter = params['ter'], z = params['z'],
        ntrials = round(ntrials*dim(obs)[1]/length(drift)), s = sigma, dt = dt, postdriftmod = params['vratio'], 
        t2distribution = rep(obs[,confRTname],times=ntrials)
      ))
      names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
      predictions['drift'] <- drift[d]
      predictions[diffname] <- difficulty[d]
    }
    else{
      temp_pred <- data.frame(DDM_with_confidence_slow_fullconfRT(
        v = drift[d], a = params['a'], ter = params['ter'], z = params['z'],
        ntrials = round(ntrials*dim(obs)[1]/length(drift)), s = sigma, dt = dt, postdriftmod = params['vratio'], 
        t2distribution = rep(obs[,confRTname],times=ntrials)
      ))
      names(temp_pred) <- c('rt','resp','cor','evidence2','rt2','cj')
      temp_pred['drift'] <- drift[d]
      temp_pred[diffname] <- difficulty[d]
      predictions <- fastmerge(predictions,temp_pred)
      
    }
  }
  
  #' Evaluate goodness-of-fit with RT quantiles for correct and errors
  c_pred <- predictions[predictions$cor == 1,]
  e_pred <- predictions[predictions$cor == 0,]
  
  c_obs <- obs[obs[,accname] %in% c(1,'cor','correct'),]
  e_obs <- obs[obs[,accname] %in% c(-1,0,'err','error'),]
  
  obs_props <- NULL; pred_props <- NULL
  for (d in 1:length(drift)) {
    #' Now, get the quantile RTs on the "obs data" for correct and error distributions 
    #' separately (for quantiles .1, .3, .5, .7, .9)
    c_quantiles <- quantile(c_obs[c_obs[,diffname] == difficulty[d],]$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    e_quantiles <- quantile(e_obs[e_obs[,diffname] == difficulty[d],]$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    if (any(is.na(e_quantiles))) {
      e_quantiles <- rep(0,5)
    }
    if (any(is.na(c_quantiles))) {
      c_quantiles <- rep(0,5)
    }      
    #' to combine correct and incorrect we scale the expected interquantile 
    #' probability by the proportion of correct and incorrect respectively
    prop_obs_c <- dim(c_obs[c_obs[,diffname] == difficulty[d],])[1] / dim(obs)[1]
    prop_obs_e <- dim(e_obs[e_obs[,diffname] == difficulty[d],])[1] / dim(obs)[1]
    
    c_obs_proportion = prop_obs_c * c(.1, .2, .2, .2, .2, .1)
    e_obs_proportion = prop_obs_e * c(.1, .2, .2, .2, .2, .1)
    obs_props <- c(obs_props,c_obs_proportion,e_obs_proportion)
    
    c_pred_rt <- sort(c_pred[c_pred[,diffname] == difficulty[d],]$rt)
    e_pred_rt <- sort(e_pred[e_pred[,diffname] == difficulty[d],]$rt)
    # now, get the proportion of responses that fall between the obs quantiles when applied to the pred data (scale by N?)
    c_pred_proportion <- c(
      sum(c_pred_rt <= c_quantiles[1]),
      sum(c_pred_rt <= c_quantiles[2]) - sum(c_pred_rt <= c_quantiles[1]),
      sum(c_pred_rt <= c_quantiles[3]) - sum(c_pred_rt <= c_quantiles[2]),
      sum(c_pred_rt <= c_quantiles[4]) - sum(c_pred_rt <= c_quantiles[3]),
      sum(c_pred_rt <= c_quantiles[5]) - sum(c_pred_rt <= c_quantiles[4]),
      sum(c_pred_rt > c_quantiles[5])
    ) / dim(predictions)[1]
    
    e_pred_proportion <- c(
      sum(e_pred_rt <= e_quantiles[1]),
      sum(e_pred_rt <= e_quantiles[2]) - sum(e_pred_rt <= e_quantiles[1]),
      sum(e_pred_rt <= e_quantiles[3]) - sum(e_pred_rt <= e_quantiles[2]),
      sum(e_pred_rt <= e_quantiles[4]) - sum(e_pred_rt <= e_quantiles[3]),
      sum(e_pred_rt <= e_quantiles[5]) - sum(e_pred_rt <= e_quantiles[4]),
      sum(e_pred_rt > e_quantiles[5])
    ) / dim(predictions)[1]
    pred_props <- c(pred_props,c_pred_proportion,e_pred_proportion)
  }
  ddm_cost <- sum((obs_props-pred_props)^2)
  
  return(ddm_cost)
}
ldc.nn.fit <- function(params,obs,ddm_params,dt=.001,sigma=0.1,
                       Nupdate_per_trial=1000,w0=c(10,-1,1),returnFit=T,
                       confRTname="RTconf",diffname="difflevel",respname="resp",
                       totRTname='rt2',targetname='cj',accname='cor',beta_input=.1,
                       error_type='mse',binning=F,nbin=6){
  #' Step 1 : Use DDM bound and drift rate to infer evidence accumulated at each trial
  #' Step 2 : Gradiant descent 
  #' Step 3 : Global cost is DDM cost + NN cost
  #' Loop and optimize NN hyperparameters (learning rate + batch size)
  
  drift <- ddm_params[2:length(ddm_params)] # Make it more flexible
  bound <- ddm_params[1]
  difficulty <- sort(unique(obs[,diffname]))
  
  #' Predict single trial amount of evidence accumulated
  obs['evidence'] <- bound
  
  if (!("trial" %in% names(obs)) ) {
    obs$trial <- 1:nrow(obs)
  }
  #' We repeat each trial several times to have the ANN update parameters
  #' several times per trial and also to take into account the stochastic nature
  #' of estimating single trial evidence accumulation process 
  obs_nn <- obs[rep(seq_len(nrow(obs)), each=Nupdate_per_trial), ]
  
  for (trial in seq(1,dim(obs_nn)[1],Nupdate_per_trial)) {
    #' Post decision drift rate sign depends on accuracy 
    if (obs_nn[trial,accname] %in% c(1,'correct','cor')) {
      obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] <- 
        obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] + 
        DDM_fixed_time(v = drift[difficulty==obs_nn[trial,diffname]],
                       time=obs_nn[trial,confRTname],ntrials=Nupdate_per_trial,s=sigma,dt=dt)[,1]
    }else if (obs_nn[trial,accname] %in% c(-1,0,'error','err')) {
      obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] <- 
        obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] + 
        DDM_fixed_time(v = - drift[difficulty==obs_nn[trial,diffname]],
                       time=obs_nn[trial,confRTname],ntrials=Nupdate_per_trial,s=sigma,dt=dt)[,1]
    }
  }
  
  #' Input data (evidence and intercept)
  x = matrix(c(obs_nn$evidence, #ev
               rep(beta_input,dim(obs_nn)[1]), #bias
               1/sqrt(obs_nn[,totRTname])),  ncol=3) #time
  
  #' Output (confidence)
  y = obs_nn[,targetname] 
  
  #' Initialize weights
  w <- w0
  
  results <- train_model(x,w,y,eta=params,error_type = error_type,trace=F,
                         binning=binning,nbin=nbin)
  
  if (returnFit) {
    return(results$err)
  }else{
    if (!("cj_pred" %in% names(obs_nn))) {
      obs_nn$cj_pred <- NA
    }
    results <- train_model(x,w,y,eta=params,error_type = error_type,trace=T,
                           binning=binning,nbin=nbin)
    trial_weight <- results$trace[seq(Nupdate_per_trial,dim(obs_nn)[1],Nupdate_per_trial),]
    for (trial in 1:dim(trial_weight)[1]) {
      obs_nn[(Nupdate_per_trial*(trial-1)+1):(Nupdate_per_trial*trial),"cj_pred"] <- 
        combine_input(x[(Nupdate_per_trial*(trial-1)+1):(Nupdate_per_trial*trial),],
                      c(trial_weight[trial,c(1,2)],1),binning=binning,nbin=nbin)
    }
    # obs_nn$cj_pred <- combine_input(x,results$w,binning=binning,nbin=nbin)
    y_pred <- with(obs_nn,aggregate(cj_pred,by=list(trial),mean))$x
    return(y_pred)
  }
}

ldc.nn.fit.w <- function(params,obs,ddm_params,dt=.001,sigma=0.1,Nsim_error=1000,
                         Nupdate_per_trial=1000,returnFit=T,estimate_evidence = T,
                         confRTname="RTconf",diffname="difflevel",respname="resp",
                         totRTname='rt2',targetname='cj',accname='cor',beta_input=.1,
                         error_type1='cross-entropy',error_type2='mse',binning=F,nbin=6,
                         shuffle=T,cost="separated",aggreg_pred="mean"){
  #' Step 1 : Use DDM bound and drift rate to infer evidence accumulated at each trial
  #' Step 2 : Gradiant descent 
  #' Step 3 : Global cost is DDM cost + NN cost
  #' Loop and optimize NN hyperparameters (learning rate + batch size)
  
  if (estimate_evidence) {
    drift <- ddm_params[2:length(ddm_params)] # Make it more flexible
    bound <- ddm_params[1]
    difficulty <- sort(unique(obs[,diffname]))
    
    #' Predict single trial amount of evidence accumulated
    obs['evidence'] <- bound
    
    if (!("trial" %in% names(obs)) ) {
      obs$trial <- 1:nrow(obs)
    }
    #' We repeat each trial several times to have the ANN update parameters
    #' several times per trial and also to take into account the stochastic nature
    #' of estimating single trial evidence accumulation process 
    obs_nn <- obs[rep(seq_len(nrow(obs)), each=Nupdate_per_trial), ]
    obs_err <- obs[rep(seq_len(nrow(obs)), each=Nsim_error), ]
    
    for (trial in seq(1,dim(obs_nn)[1],Nupdate_per_trial)) {
      #' Post decision drift rate sign depends on accuracy 
      if (obs_nn[trial,accname] %in% c(1,'correct','cor')) {
        obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] <- 
          obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] + 
          DDM_fixed_time(v = drift[difficulty==obs_nn[trial,diffname]],
                         time=obs_nn[trial,confRTname],ntrials=Nupdate_per_trial,s=sigma,dt=dt)[,1]
      }else if (obs_nn[trial,accname] %in% c(-1,0,'error','err')) {
        obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] <- 
          obs_nn[trial:(trial+Nupdate_per_trial-1),'evidence'] + 
          DDM_fixed_time(v = - drift[difficulty==obs_nn[trial,diffname]],
                         time=obs_nn[trial,confRTname],ntrials=Nupdate_per_trial,s=sigma,dt=dt)[,1]
      }
    }
    for (trial in seq(1,dim(obs_err)[1],Nsim_error)) {
      #' Post decision drift rate sign depends on accuracy 
      if (obs_err[trial,accname] %in% c(1,'correct','cor')) {
        obs_err[trial:(trial+Nsim_error-1),'evidence'] <- 
          obs_err[trial:(trial+Nsim_error-1),'evidence'] + 
          DDM_fixed_time(v = drift[difficulty==obs_err[trial,diffname]],
                         time=obs_err[trial,confRTname],ntrials=Nsim_error,s=sigma,dt=dt)[,1]
      }else if (obs_err[trial,accname] %in% c(-1,0,'error','err')) {
        obs_err[trial:(trial+Nsim_error-1),'evidence'] <- 
          obs_err[trial:(trial+Nsim_error-1),'evidence'] + 
          DDM_fixed_time(v = - drift[difficulty==obs_err[trial,diffname]],
                         time=obs_err[trial,confRTname],ntrials=Nsim_error,s=sigma,dt=dt)[,1]
      }
    }
    
  }else{
    obs_nn <- obs
    if (shuffle) {
      obs_nn <- do.call(rbind, 
                        lapply(split(obs_nn, obs_nn$trial), 
                               function(x) x[sample(nrow(x), Nupdate_per_trial), ]))
    }
  }
  
  #' Input data (evidence and intercept)
  x = matrix(c(obs_nn$evidence, #ev
               rep(beta_input,dim(obs_nn)[1]), #bias
               1/sqrt(obs_nn[,totRTname])),  ncol=3) #time

  x_err = matrix(c(obs_err$evidence, #ev
               rep(beta_input,dim(obs_err)[1]), #bias
               1/sqrt(obs_err[,totRTname])),  ncol=3) #time
  
  #' Output (confidence)
  y = obs_nn[,targetname] 
  
  #' Initialize weights
  w <- params[1:3]
  
  results <- train_model(x,w,y,eta=params[4],error_type1 = error_type1,trace=F,
                         binning=binning,nbin=nbin,Nupdate_per_trial = Nupdate_per_trial,
                         cost=cost,x_err = x_err,Nsim_error=Nsim_error,
                         error_type2 = error_type2)
  
  if (returnFit) {
    return(results$err)
  }else{
    if (!("cj_pred" %in% names(obs_nn))) {
      obs_nn$cj_pred <- NA
    }
    results <- train_model(x,w,y,eta=params[4],error_type1 = error_type1,trace=T,
                           Nupdate_per_trial = Nupdate_per_trial,binning=binning,
                           nbin=nbin,cost=cost,error_type2 = error_type2,x_err=x_err)
    
    trial_weight <- results$trace[seq(Nupdate_per_trial,dim(obs_nn)[1],Nupdate_per_trial),]
    for (trial in 1:dim(trial_weight)[1]) {
      obs_nn[(Nupdate_per_trial*(trial-1)+1):(Nupdate_per_trial*trial),"cj_pred"] <- 
        combine_input(matrix(x[(Nupdate_per_trial*(trial-1)+1):(Nupdate_per_trial*trial),],ncol=3),
                      c(trial_weight[trial,c(1,2)],1),binning=binning,nbin=nbin)
    }
    # obs_nn$cj_pred <- combine_input(x,results$w,binning=binning,nbin=nbin)
    if (aggreg_pred=="mean") {
      y_pred <- with(obs_nn,aggregate(cj_pred,by=list(trial),mean))$x
    } else if (aggreg_pred=="mode"){
      y_pred <- with(obs_nn,aggregate(cj_pred,by=list(trial),Mode))$x
    }
    return(list(pred=y_pred,trace=trial_weight))
  }
}

ldc.fit <- function(params,ddm_params1,ddm_params2,obs1,obs2,dt=.001,sigma=0.1,
                    Nupdate_per_trial=1000,returnFit=T,Nsim_error=1000,
                    confRTname="RTconf",diffname="difflevel",respname="resp",
                    totRTname='rt2',targetname='cj',accname='cor',beta_input=.1,
                    error_type1='cross-entropy',error_type2='mse',binning=F,nbin=6,
                    cost='separated'){
  fit1 <- ldc.nn.fit.w(params[1:4],obs1,ddm_params1,dt=dt,sigma=sigma,
                       Nupdate_per_trial=Nupdate_per_trial,returnFit=returnFit,
                       confRTname=confRTname,diffname=diffname,respname=respname,
                       totRTname=totRTname,targetname=targetname,accname=accname,
                       beta_input=beta_input,error_type1=error_type1,binning=binning,nbin=nbin,
                       Nsim_error=Nsim_error,cost=cost,error_type2=error_type2)
  fit2 <- ldc.nn.fit.w(params[5:8],obs2,ddm_params2,dt=dt,sigma=sigma,
                       Nupdate_per_trial=Nupdate_per_trial,returnFit=returnFit,
                       confRTname=confRTname,diffname=diffname,respname=respname,
                       totRTname=totRTname,targetname=targetname,accname=accname,
                       beta_input=beta_input,error_type1=error_type1,binning=binning,nbin=nbin,
                       Nsim_err=Nsim_error,cost=cost,error_type2=error_type2)
  if (returnFit) {
    return(fit1 + fit2)
  }else{
    return(c(fit1,fit2))
  }
}
