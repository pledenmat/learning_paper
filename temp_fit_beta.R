train_beta <- read.csv("train_beta.csv")
train_beta$rt2 <- train_beta$rt + train_beta$RTconf

subs <- unique(train_beta$sub)
conditions <- sort(unique(train_beta$condition))
difflevels <- sort(unique(train_beta$difflevel)) # Important to sort to match with drift order


#' FIT BETA EXPERIMENT
temp_dat <- subset(train_beta,sub==s&condition==cond)
ddm_file <- paste0('fit/beta/ddmfit_',s,'_',cond,'.Rdata')
if(file.exists(ddm_file)){
  load(ddm_file)
}else{
  optimal_params <- DEoptim(ddm.fit, # function to optimize
                            lower = params_lower, upper = params_upper,obs = temp_dat,
                            dt = dt, sigma = s,ntrials=20,
                            control=c(itermax=1000,steptol=100,reltol=.001,NP=50))
  ddm.results <- summary(optimal_params)
  #save individual results
  save(ddm.results, file=ddm_file)
}
ddm_params <- ddm.results$optim$bestmem[c(1,5:length(ddm.results$optim$bestmem))]
ldc_file <- paste0('fit/beta/ldcfit_',s,'_',cond,'.Rdata')
if (file.exists(ldc_file)) {
  load(ldc_file)
}else{
  optimal_params <- DEoptim(ldc.nn.fit,ddm_params=ddm_params, obs = temp_dat,
                            lower = 0, upper = 1,
                            dt = dt, sigma = sigma,
                            control=c(itermax=1000,steptol=100,reltol=.001,NP=10))
  ldc.results <- summary(optimal_params)
  #save individual results
  save(ldc.results, file=ldc_file)
}
par_b[par_b$sub==s&par_b$condition==cond,"bound"] <- ddm.results$optim$bestmem[1]
par_b[par_b$sub==s&par_b$condition==cond,"ter"] <- ddm.results$optim$bestmem[2]
par_b[par_b$sub==s&par_b$condition==cond&par_b$difflevel==difflevels[1],"drift"] <- ddm.results$optim$bestmem[5]
par_b[par_b$sub==s&par_b$condition==cond&par_b$difflevel==difflevels[2],"drift"] <- ddm.results$optim$bestmem[6]
par_b[par_b$sub==s&par_b$condition==cond&par_b$difflevel==difflevels[3],"drift"] <- ddm.results$optim$bestmem[7]
par_b[par_b$sub==s&par_b$condition==cond,"eta"] <- ldc.results$optim$bestmem[1]
par_b[par_b$sub==s&par_b$condition==cond,"cost_ddm"] <- ddm.results$optim$bestval
par_b[par_b$sub==s&par_b$condition==cond,"cost_ldc"] <- ldc.results$optim$bestval