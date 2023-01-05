rm(list=ls())
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
library(myPackage)
source("rw_createHM_driftsign.r")
library(reshape)
library(timeSeries)
library(fields)
library(MALDIquant)

model_conf <- function(a,b,data){
  return(( (1+data$resp)/2)*(1 /(1 + exp( (1/sqrt(data$rt2+.00001))* (-a*data$evidence2 - b))))  + ((1-data$resp)/2)*(1 /(1 + exp((1/sqrt(data$rt2+.00001))* (a*data$evidence2 - b))) ))
}

# For plotting ------------------------------------------------------------


error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
col_low <- rgb(237,248,177,maxColorValue = 255)
col_med <- rgb(127,205,187,maxColorValue = 255)
col_high <- rgb(44,127,184,maxColorValue = 255)
lwd_trace <- 3

# Parameters --------------------------------------------------------------
v <- .1 #' Drift rate
s <- .1 #' Within-trial noise
a <- .07 #' Boundary distance from 0
ter <- .4 #' Non-decision time
z <- 0 #' Starting point bias
dt <- .001 #' Time step
t2time <- .5 #' Post-decision accumulation time
vratio <- 1
ntrials <- 100000
params <- c(a,ter,z,ntrials,s,dt,t2time,vratio,v)

ev_bound <- .5; ev_window <- dt*10; upperRT <- 5
ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
timesteps <- upperRT/dt
mu <- c(-v,v)
nsim <- 100000/length(mu)

# Data generation ---------------------------------------------------------
# We'll generate some data that we can use to train the perceptron
Data <- chi_square_optim_DDM(params,returnFit = 0,observations = NULL)
Data$evidence2 <- Data$evidence2*Data$resp
Data$resp <- 1
# Train function --------------------------------------------------------------
model_linear <- function(x,w){
  y_j <- as.numeric(x %*% w)
  pred_j <- 1/(1+exp(-y_j))
  return(pred_j)
}
model_nolinear <- function(x,w){
  if(is.vector(x)){
    y_j <- as.numeric(x[3]*w[3]*(x[c(1,2)] %*% w[c(1,2)]))
  }else{
    y_j <- as.numeric(x[,3]*w[3]*(x[,c(1,2)] %*% w[c(1,2)]))
  }
  pred_j <- 1/(1+exp(-y_j))
  return(pred_j)
}
cross_entropy <- function(y,y_pred){
  return(- sum( (1 - y)*log(1 - y_pred) + y*log(y_pred))/length(y))
}
grad <- function(x,y,y_pred){
  return((y - y_pred)*x[,c(1,2)]/sqrt(x[,3]))
}
train_model <- function(x,w,y,eta=.1,momentum=0,
                        M = 1,th=0.0000001,verbose=F,trace=F,model="nolinear"){
  #' Add test that x,w and y dimensions are compatible
  ntrain <- length(y)
  dat_iter <- data.frame(alpha = rep(NA,M*ntrain),
                         beta = rep(NA,M*ntrain),
                         err = rep(NA,M*ntrain)) 
  for (i in 1:M){
    # print(paste('Epoch starts: ', i))
    
    ## We reshuffle the order of the datapoint for each epoch.
    index = 1:ntrain
    # index = sample(index)
    
    for (j in index){
      
      if(model=="linear"){
        pred_j <- model_linear(x[j,],w)
      }else if (model == "nolinear") {
        pred_j <- model_nolinear(x[j,],w)
      }
      # Update weights (time is fixed)
      w[c(1,2)] = w[c(1,2)] + eta*grad(matrix(x[j,],ncol=3),y[j],pred_j)
      
      if (verbose == T){
        print(paste('  -> updating data point ', j, ' : '))
        print(paste('     -> alpha: ' ,w[1]))
        print(paste('     -> beta: ' ,w[2]))
      }
      if(model=="linear"){
        y_pred <- model_linear(x,w)
      }else if (model == "nolinear") {
        y_pred <- model_nolinear(x,w)
        y_pred[y_pred==1] <- .9999999 # Avoid log(0)
        y_pred[y_pred==0] <- .0000001
      }
      err = cross_entropy(y,y_pred)
      if (trace) {
        dat_iter[(i-1)*ntrain+j,"alpha"] <- w[1]
        dat_iter[(i-1)*ntrain+j,"beta"] <- w[2]
        dat_iter[(i-1)*ntrain+j,"err"] <- err
      }
    }  
    # print(paste('Epoch ends: ', i, ' WITH accuracy: ', acc))
    if (err <= th){
      # break
    }
  }
  results <- list(w=w,err=err,param_iter=dat_iter,niter=i*ntrain)
  return(results)
}

# Train with trace --------------------------------------------------------
#' Input data (evidence and intercept)
x = matrix(c(Data$evidence2, #ev
             rep(1,dim(Data)[1]), #bias
             1/sqrt(Data$rt2)),  ncol=3) #time

#' Output (confidence)
y = Data$cor 

#' Initialize weights
w <- c(0,-1,1)

M <- 1
eta <- .1
if (!file.exists("train_trace_loweta.Rdata")) {
  results_loweta <- train_model(x,w,y,eta=eta,M=M,trace=T,verbose=T)
  save(results_loweta,file="train_trace_loweta.Rdata")
}else{
  load("train_trace_loweta.Rdata")
}

#' Try other learning rates
eta <- .9
if (!file.exists("train_trace_higheta.Rdata")) {
  results_higheta <- train_model(x,w,y,eta=eta,M=M,trace=T,verbose=T)
  save(results_higheta,file="train_trace_higheta.Rdata")
}else{
  load("train_trace_higheta.Rdata")
}

eta <- .5
if (!file.exists("train_trace_medeta.Rdata")) {
  results_medeta <- train_model(x,w,y,eta=eta,M=M,trace=T,verbose=T)
  save(results_medeta,file="train_trace_medeta.Rdata")
}else{
  load("train_trace_medeta.Rdata")
}


#' Initialize with higher alpha
w <- c(50,-1,1)

eta <- .1
if (!file.exists("train_trace_loweta_highalpha.Rdata")) {
  results_loweta_highalpha <- train_model(x,w,y,eta=eta,M=M,trace=T,verbose=T)
  save(results_loweta_highalpha,file="train_trace_loweta_highalpha.Rdata")
}else{
  load("train_trace_loweta_highalpha.Rdata")
}


eta <- .5
if (!file.exists("train_trace_medeta_highalpha.Rdata")) {
  results_medeta_highalpha <- train_model(x,w,y,eta=eta,M=M,trace=T,verbose=T)
  save(results_medeta_highalpha,file="train_trace_medeta_highalpha.Rdata")
}else{
  load("train_trace_medeta_highalpha.Rdata")
}

eta <- .9
if (!file.exists("train_trace_higheta_highalpha.Rdata")) {
  results_higheta_highalpha <- train_model(x,w,y,eta=eta,M=M,trace=T,verbose=T)
  save(results_higheta_highalpha,file="train_trace_higheta_highalpha.Rdata")
}else{
  load("train_trace_higheta_highalpha.Rdata")
}


#' Plot the trace
plot(results_loweta$param_iter$alpha,type='l',col=col_low,lwd=lwd_trace,
     ylim=c(0,50),xlab="Time (trials)",ylab="Alpha",frame=F)
lines(results_higheta$param_iter$alpha,col=col_high,lwd=lwd_trace)
lines(results_medeta$param_iter$alpha,col=col_med,lwd=lwd_trace)
lines(results_medeta_highalpha$param_iter$alpha,col=col_med,lwd=lwd_trace)
lines(results_loweta_highalpha$param_iter$alpha,col=col_low,lwd=lwd_trace)
lines(results_higheta_highalpha$param_iter$alpha,col=col_high,lwd=lwd_trace)
legend("topright",legend=c("0.1","0.5","0.9"),col=c(col_low,col_med,col_high),
       lty=1,title="Learning rate",bty='n',lwd=lwd_trace)

plot(results_loweta$param_iter$beta,type='l',col=col_low,lwd=lwd_trace,
     xlab="Time (trials)",ylab="Beta",ylim=c(-3.5,3.5))
lines(results_higheta$param_iter$beta,col=col_high,lwd=lwd_trace)
lines(results_medeta$param_iter$beta,col=col_med,lwd=lwd_trace)
lines(results_medeta_highalpha$param_iter$beta,col=col_med,lwd=lwd_trace)
lines(results_loweta_highalpha$param_iter$beta,col=col_med,lwd=lwd_trace)
lines(results_higheta_highalpha$param_iter$beta,col=col_med,lwd=lwd_trace)

plot(results_loweta$param_iter$err,type='l',col=col_low,lwd=lwd_trace,
     xlab="Time (trials)",ylab="Cross-entropy")
lines(results_higheta$param_iter$err,col=col_high,lwd=lwd_trace)
lines(results_medeta$param_iter$err,col=col_med,lwd=lwd_trace)
lines(results_medeta_highalpha$param_iter$err,col=col_med,lwd=lwd_trace)
lines(results_loweta_highalpha$param_iter$err,col=col_med,lwd=lwd_trace)
lines(results_higheta_highalpha$param_iter$err,col=col_med,lwd=lwd_trace)


#' Let's try with beta input of .1 instead of 1
x = matrix(c(Data$evidence2, #ev
             rep(.1,dim(Data)[1]), #bias
             1/sqrt(Data$rt2)),  ncol=3) #time
w <- c(0,-1,1)
eta <- .1
if (!file.exists("train_trace_loweta_betainput.Rdata")) {
  results_loweta_betainput <- train_model(x,w,y,eta=eta,M=M,trace=T,verbose=T)
  save(results_loweta_betainput,file="train_trace_loweta_betainput.Rdata")
}else{
  load("train_trace_loweta_betainput.Rdata")
}

plot(results_loweta_betainput$param_iter$alpha,type='l',lwd=lwd_trace,
     ylim=c(0,50),xlab="Time (trials)",ylab="Alpha",frame=F,main = "Beta input set to 0.1")
plot(results_loweta_betainput$param_iter$beta,type='l',lwd=lwd_trace,frame=F,
     xlab="Time (trials)",ylab="Beta")
plot(results_loweta_betainput$param_iter$err,type='l',lwd=lwd_trace,frame=F,
     xlab="Time (trials)",ylab="Cross-entropy")
# Test starting weights -----------------------------------------------------------
n_train <- ntrials
temp_dat <- Data[sample(1:n_train),]
#' Input data (evidence and intercept)
x = matrix(c(temp_dat$evidence2, #ev
             rep(.1,dim(temp_dat)[1]), #bias
             1/sqrt(temp_dat$rt2)),  ncol=3) #time

#' Output (confidence)
y = temp_dat$cor 


#' Initialize weights
n_start <- 10
a0 <- rep(seq(0,50,length.out=n_start),each=n_start)
b0 <- rep(seq(-2,2,length.out=n_start),n_start)
w = matrix(c(a0,b0,rep(1,n_start^2)),ncol=3)
# c(a0,b0,1)

M = 1000000/n_train        # number of epochs to run
eta = 0.1       # learning rate
th = 0.00000001          # threshold to stop
verbose = F   # whether detailed weight update info is printed
model = "nolinear"
momentum <- 0

results_w0 <- data.frame(alpha=rep(NA,dim(w)[1]),
                         beta=rep(NA,dim(w)[1]),
                         CE=rep(NA,dim(w)[1]),
                         niter=rep(NA,dim(w)[1]))
for (iter in 1:dim(w)[1]) {
  print(paste("Training model",iter,"of",dim(w)[1]))
  results <- train_model(x,w[iter,],y,M=M,th=th)
  results_w0[iter,"alpha"] <- results$w[1]
  results_w0[iter,"beta"] <- results$w[2]
  results_w0[iter,"CE"] <- results$acc
  results_w0[iter,"niter"] <- results$niter
}

#Threshold was always reached
plot(results_w0$CE,ylab="CE",ylim=c(min(min(results_w0$CE,th)),max(max(results_w0$CE),th)),frame=F)
abline(h=th,col="red")

#Plot number of iteration before reaching threshold
plot(results_w0$niter/500,frame=F,xaxt="n",xlab="Alpha",ylab="N epochs")
axis(1,at=seq(0,100,length.out=10),labels = round(seq(0,50,length.out=n_start)))
abline(v=alpha*2)
print(w)

plot(results_w0$alpha,frame=F,ylab="Alpha")
abline(h=alpha,col="red")
plot(results_w0$beta,frame=F,ylab="Beta")
abline(h=beta,col="red")


# Test number of train trials ---------------------------------------------

n_train <- c(5,10,20,50,100,200,500)
n_rep <- 10

results_n_train <- data.frame(alpha=rep(NA,length(n_train)*n_rep),
                              beta=rep(NA,length(n_train)*n_rep),
                              CE=rep(NA,length(n_train)*n_rep),
                              niter=rep(NA,length(n_train)*n_rep),
                              n_train=rep(n_train,each=n_rep),
                              rep=rep(1:n_rep,length(n_train)))
th <- 0.00000001

for (iter in 1:length(n_train)) {
  
  M = 10000000/n_train[iter]        # number of epochs to run
  print(paste("Training model",iter,"of",length(n_train)))
  for (rep in 1:n_rep) {
    temp_dat <- Data[sample(1:dim(Data)[1],size=n_train[iter]),]
    #' Input data (evidence, bias and time)
    x = matrix(c(temp_dat$evidence2, #ev
                 rep(.1,dim(temp_dat)[1]), #bias
                 1/sqrt(temp_dat$rt2)),  ncol=3) #time
    
    #' Output (confidence)
    y = temp_dat$cj 
    
    
    #' Initialize weights
    w <- c(0,0,1)
    
    print(paste("repetition",rep,"of",n_rep))
    results <- train_model(x,w,y,M=M,th=th)
    results_n_train[(iter-1)*n_rep+rep,"alpha"] <- results$w[1]
    results_n_train[(iter-1)*n_rep+rep,"beta"] <- results$w[2]
    results_n_train[(iter-1)*n_rep+rep,"CE"] <- results$acc
    results_n_train[(iter-1)*n_rep+rep,"niter"] <- results$niter
  }
}
ce <- cast(results_n_train,rep~n_train,value = "CE")
#Threshold was always reached
plot(colMeans(ce),ylab="CE",type='p',pch=16,
     ylim=c(min(min(results_n_train$CE,th)),max(max(results_n_train$CE),th)),frame=F)
error.bar(1:length(n_train),colMeans(ce),colSds(ce,na.rm=T),lwd=3,length=0)
abline(h=th,col="red")

#Plot number of iteration before reaching threshold
niter <- cast(results_n_train,rep~n_train,value = "niter")
plot(colMeans(niter),frame=F,xaxt="n",xlab="Ntrials",ylab="Nupdate",pch=16)
error.bar(1:length(n_train),colMeans(niter),colSds(niter,na.rm=T),lwd=3,length=0)
axis(1,at=1:length(n_train),labels = n_train)

alphas <- cast(results_n_train,rep~n_train,value = "alpha")
plot(colMeans(alphas),frame=F,ylab="Alpha",pch=16,xaxt='n',xlab="Ntrials",
     ylim=c(min(min(results_n_train$alpha,alpha)),max(max(results_n_train$alpha),alpha)))
error.bar(1:length(n_train),colMeans(alphas),colSds(alphas,na.rm=T),lwd=3,length=0)
axis(1,at=1:length(n_train),labels = n_train)
abline(h=alpha,col="red")

betas <- cast(results_n_train,rep~n_train,value = "beta")
plot(colMeans(betas)/10,frame=F,ylab="Beta",pch=16,xaxt='n',xlab="Ntrials",
     ylim=c(min(min(results_n_train$beta,beta)),max(max(results_n_train$beta),beta)))
error.bar(1:length(n_train),colMeans(betas)/10,colSds(betas,na.rm=T),lwd=3,length=0)
axis(1,at=1:length(n_train),labels = n_train)
abline(h=beta,col="red")


# Test gradient map alpha/beta  -----------------------------------
alpha_min <- 0
alpha_max <- 50
alpha_step <- .5
beta_min <- -3.5
beta_max <- 3.5
beta_step <- .1
alpha_grid <- seq(alpha_min,alpha_max,alpha_step)
beta_grid <- seq(beta_min,beta_max,beta_step)
n_sample <- 1000
temp_dat <- Data[sample(1:dim(Data)[1],n_sample),]
x = matrix(c(temp_dat$evidence2, #ev
             rep(1,dim(temp_dat)[1]), #bias
             1/sqrt(temp_dat$rt2)),  ncol=3) #time
y = temp_dat$cor

grid <- data.frame(alpha=rep(alpha_grid,each=length(beta_grid)),
                   beta=rep(beta_grid,length(alpha_grid)),
                   ce=NA)
for (a in alpha_grid) {
  for (b in beta_grid) {
    w <- c(a,b,1)
    y_pred <- model_nolinear(x,w)
    grid[grid$alpha==a&grid$beta==b,"ce"] <- cross_entropy(y,y_pred)
    
  }
}
grid_mat <- as.matrix(cast(grid,alpha~beta))
grid_vec <- as.vector(grid_mat)
grid_min <- which(grid_mat == min(grid_mat), arr.ind = TRUE)


# grid_mat <- cast(grid,alpha~beta)
results_higheta$param_iter$a_grid <- match.closest(results_higheta$param_iter$alpha,alpha_grid)
results_higheta$param_iter$b_grid <- match.closest(results_higheta$param_iter$beta,beta_grid)
results_medeta$param_iter$a_grid <- match.closest(results_medeta$param_iter$alpha,alpha_grid)
results_medeta$param_iter$b_grid <- match.closest(results_medeta$param_iter$beta,beta_grid)
results_loweta$param_iter$a_grid <- match.closest(results_loweta$param_iter$alpha,alpha_grid)
results_loweta$param_iter$b_grid <- match.closest(results_loweta$param_iter$beta,beta_grid)
#' Plot with big learning rate
image.plot(1:dim(grid_mat)[1],1:dim(grid_mat)[2],grid_mat,ylab="Beta",
           xlab="Alpha",legend.shrink = .5, axes = F,main="eta = 0.9")
axis(1,at=seq(1,length(alpha_grid),round(length(alpha_grid)/10)),labels=seq(alpha_min,alpha_max,(alpha_max-alpha_min)/10));
axis(2,at=seq(1,length(beta_grid),round(length(beta_grid)/10)),labels=seq(beta_min,beta_max,(beta_max-beta_min)/10));
points(type='l',x=results_higheta$param_iter$a_grid,
       y=results_higheta$param_iter$b_grid)
points(type='p',col="red",x=grid_min[1],y=grid_min[2])


#' Plot with medium learning rate
image.plot(1:dim(grid_mat)[1],1:dim(grid_mat)[2],grid_mat,ylab="Beta",
           xlab="Alpha",legend.shrink = .5, axes = F,main="eta = 0.5")
axis(1,at=seq(1,length(alpha_grid),round(length(alpha_grid)/10)),labels=seq(alpha_min,alpha_max,(alpha_max-alpha_min)/10));
axis(2,at=seq(1,length(beta_grid),round(length(beta_grid)/10)),labels=seq(beta_min,beta_max,(beta_max-beta_min)/10));
points(type='l',x=results_medeta$param_iter$a_grid,
       y=results_medeta$param_iter$b_grid)
points(type='p',col="red",x=grid_min[1],y=grid_min[2])

#' Plot with low learning rate
image.plot(1:dim(grid_mat)[1],1:dim(grid_mat)[2],grid_mat,ylab="Beta",
           xlab="Alpha",legend.shrink = .5, axes = F,main="eta = 0.1")
axis(1,at=seq(1,length(alpha_grid),round(length(alpha_grid)/10)),labels=seq(alpha_min,alpha_max,(alpha_max-alpha_min)/10));
axis(2,at=seq(1,length(beta_grid),round(length(beta_grid)/10)),labels=seq(beta_min,beta_max,(beta_max-beta_min)/10));
points(type='l',x=results_loweta$param_iter$a_grid,
       y=results_loweta$param_iter$b_grid)
points(type='p',col="red",x=grid_min[1],y=grid_min[2])

# Grid for beta input of .1 -----------------------------------------------

alpha_min <- 0
alpha_max <- 30
alpha_step <- .2
beta_min <- -1.5
beta_max <- 6.5
beta_step <- .1
alpha_grid <- seq(alpha_min,alpha_max,alpha_step)
beta_grid <- seq(beta_min,beta_max,beta_step)
n_sample <- 1000
temp_dat <- Data[sample(1:dim(Data)[1],n_sample),]
x = matrix(c(temp_dat$evidence2, #ev
             rep(.1,dim(temp_dat)[1]), #bias
             1/sqrt(temp_dat$rt2)),  ncol=3) #time
y = temp_dat$cor

grid <- data.frame(alpha=rep(alpha_grid,each=length(beta_grid)),
                   beta=rep(beta_grid,length(alpha_grid)),
                   ce=NA)
for (a in alpha_grid) {
  for (b in beta_grid) {
    w <- c(a,b,1)
    y_pred <- model_nolinear(x,w)
    grid[grid$alpha==a&grid$beta==b,"ce"] <- cross_entropy(y,y_pred)
    
  }
}
grid_mat <- as.matrix(cast(grid,alpha~beta))
grid_vec <- as.vector(grid_mat)
grid_min <- which(grid_mat == min(grid_mat), arr.ind = TRUE)
# grid_mat <- cast(grid,alpha~beta)
results_loweta_betainput$param_iter$a_grid <- match.closest(results_loweta_betainput$param_iter$alpha,alpha_grid)
results_loweta_betainput$param_iter$b_grid <- match.closest(results_loweta_betainput$param_iter$beta,beta_grid)

image.plot(1:dim(grid_mat)[1],1:dim(grid_mat)[2],grid_mat,ylab="Beta",xlab="Alpha",legend.shrink = .5, axes = F)
axis(1,at=seq(1,length(alpha_grid),round(length(alpha_grid)/10)),labels=seq(alpha_min,alpha_max,(alpha_max-alpha_min)/10));
axis(2,at=seq(1,length(beta_grid),round(length(beta_grid)/10)),labels=seq(beta_min,beta_max,(beta_max-beta_min)/10));
points(type='p',col="red",x=grid_min[1],y=grid_min[2])
points(type='l',x=results_loweta_betainput$param_iter$a_grid,
       y=results_loweta_betainput$param_iter$b_grid)
# Test with empirical data ------------------------------------------------
emp_data <- read.csv("data_fbmod.csv")
subs <- sort(unique(emp_data$sub))
condlab <- sort(unique(emp_data$condition))
n_sub <- length(subs)
n_cond <- length(condlab)
for (s in 1:n_sub) {
  for (cond in 1:n_cond) {
    temp_dat <- subset(emp_data,sub==subs[s]&condition==condlab[cond])
    
    x = matrix(c(temp_dat$evidence2, #ev
                 rep(.1,dim(temp_dat)[1]), #bias
                 1/sqrt(temp_dat$rt2)),  ncol=3) #time
    
    #' Output (confidence)
    y = temp_dat$cj 
    
  }
}

trace <- trace[complete.cases(trace),]
# plot(trace$alpha)
# plot(trace$beta)

v <- seq(0,1,.01)
w <- v - v^2
plot(w)

# Trash -------------------------------------------------------------------
#' 
#' #' Build p(cor) heatmap from drift rate used to generate data 
#' output <- RW_createHM_driftsign(mu, dt=dt, nsim=nsim, ev_bound=ev_bound, 
#'                                 ev_window=ev_window, upperRT=upperRT)
#' 
#' #' Reshape the heatmap to fit alpha & beta to it
#' hm_up <- output$upper
#' 
#' hmdf <- data.frame(cj=c(hm_up),ev=rep(ev_mapping,each=timesteps),
#'                    t=rep(seq(0,upperRT-dt,dt),length(ev_mapping)))
#' hmdf$resp <- 1
#' 
#' hmdf <- hmdf[complete.cases(hmdf),]
#' 
#' # Fit alpha and beta to heatmap
#' fit <- nls(cj ~ (1/(1+exp((1/sqrt(t+.00001)) *-(a *ev + b)))), 
#'            data=hmdf, start=list(a=15, b=0),control = nls.control(warnOnly = T),trace=T)
#' alpha <- coef(fit)[[1]] 
#' beta <- coef(fit)[[2]] 
#' 
#' #' Confidence as generated by alpha and beta
#' Data$cj <- model_conf(alpha,beta,Data)
