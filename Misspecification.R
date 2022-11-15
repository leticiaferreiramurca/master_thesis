## DIC
dic <- function(x, y, logLikelihood0, postSamples){
  
  logLikelihood <- function(theta) logLikelihood0(theta, y, x)
  
  thetaBayes <- colMeans(postSamples)
  
  dbar <- mean(-2 * apply(postSamples, 1, logLikelihood))
  dhat <- -2 * logLikelihood(thetaBayes)
  dic <- 2 * dbar - dhat
  
  return(dic)
}

## EAIC
eaic <- function(x,y, logLikelihood0, postSamples){
  
  logLikelihood <- function(theta) logLikelihood0(theta, y, x)
  
  dbar <- mean(-2 * apply(postSamples, 1, logLikelihood) )
  eaic <- dbar + 2*ncol(postSamples)
  
  return(eaic)
}

## EBIC
ebic <- function(x,y, logLikelihood0, postSamples){
  n <- length(y)
  
  logLikelihood <- function(theta) logLikelihood0(theta, y, x)
  
  dbar <- mean(-2 * apply(postSamples, 1, logLikelihood) )
  ebic <- dbar + ncol(postSamples)*log(n)
  
  return(ebic)
}


################################################################################
###### Misspecification ########################################################
################################################################################
library(rstan)
options(mc.cores = 7)
library("loo")

################################################################################
#### Estimação Double Lomax ####################################################
#################################################################################

#generate fake data
F_cauchit <- function(x, lambda){
  
  prob = pcauchy(x)^lambda
  
  return(prob)
}

n = 5000
r=100

beta0 = 0
beta1 = 1
lambda = 4

### 100 amostras
XX = matrix(ncol = r, nrow = n)
yy = matrix(ncol = r, nrow = n)
for(i in 1:r){
  set.seed(i)
  XX[,i] = runif(n, -3, 3)
  eta = beta0 + beta1*x
  p <- F_cauchit(eta, lambda) 
  yy[,i] = rbinom(n,1,p)
}

mean(apply(yy,2,mean))

########################################################
#### Estimação Logística ###########################
#########################################################
loglike_log <- function(par,y,x) {
  xbeta=x%*%par	
  p=exp(xbeta)/(1+exp(xbeta)) 	
  logvero=sum(y*log(p) + (1-y)*log(1-p))
  return(logvero)
}


model_logistic <- stan_model('dlogistic.stan')
beta0 <- vector()
beta1 <- vector()

measures <- data.frame(matrix(ncol = 5, nrow = r))
names(measures) <- c("LOO", "WAIC", "DIC", "EAIC", "EBIC")

for(i in 1:r){
  set.seed(i)
  #400 samples
  fit_logistic <- sampling(model_logistic, list(n = n, y=yy[,i], x=XX[,i]), iter = 200, chains = 4)
  
  ## Looic e Waic
  log_lik_1 <- extract_log_lik(fit_logistic)
  measures[i,1] <- loo(log_lik_1)$looic
  measures[i,2] <- waic(log_lik_1)$waic
  
  #DIC, EAIC e EBIC
  postSamples <- Reduce(cbind, extract(fit_logistic, pars=c("beta0", "beta1")))
  measures[i,3] <- dic(model.matrix(~XX[,i]), yy[,i], loglike_log, postSamples)
  measures[i,4] <- eaic(model.matrix(~XX[,i]), yy[,i], loglike_log, postSamples)
  measures[i,5] <- ebic(model.matrix(~XX[,i]), yy[,i], loglike_log, postSamples)
  print(i)
}


write.csv(measures, "misspecification_logistic25.csv")



########################################################
#### Estimação Double Lomax ###########################
#########################################################
F_dlomax <- function(x){
  prob = vector(length = length(x))
  
  prob[which(x>0)] = 1 - 0.5*(1/(1+x[x>0]))
  prob[which(x<=0)] = 0.5*(1/(1 - x[x<=0]))
  
  return(prob)
}


loglike_dlomax <- function(par,y,x){
  eta <-x %*% par
  p <- F_dlomax(eta)
  logvero=sum(y*log(p) + (1-y)*log(1-p))
  return(logvero) 
}


model_lomax <- stan_model('dlomax.stan')
beta0_lomax <- vector()
beta1_lomax <- vector()

for(i in 1:r){
  set.seed(i)
  #400 samples
  fit_dlomax <- sampling(model_lomax, list(n = n, y=yy[,i], x=XX[,i]), iter = 200, chains = 4)
  
  ## Looic e Waic
  log_lik_1 <- extract_log_lik(fit_dlomax)
  measures[i,1] <- loo(log_lik_1)$looic
  measures[i,2] <- waic(log_lik_1)$waic
  
  #DIC, EAIC e EBIC
  postSamples <- Reduce(cbind, extract(fit_dlomax, pars=c("beta0", "beta1")))
  measures[i,3] <- dic(model.matrix(~XX[,i]), yy[,i], loglike_dlomax, postSamples)
  measures[i,4] <- eaic(model.matrix(~XX[,i]), yy[,i], loglike_dlomax, postSamples)
  measures[i,5] <- ebic(model.matrix(~XX[,i]), yy[,i], loglike_dlomax, postSamples)
  print(i)
}

write.csv(measures, "misspecification_dlomax25.csv")


##############################################################################
### Estimação Potência Lomax ######################################x##############
###################################################################################

Fplomax <- function(x, loglambda){
  lambda = exp(loglambda)
  p = F_dlomax(x)^lambda
  return(p)
}

loglike_dplomax <- function(par,y,x){
  eta <-x %*% par[-length(par)]
  p <- Fplomax(eta, par[length(par)])
  logvero = sum(y*log(p) + (1-y)*log(1-p))
  return(logvero) 
}

model_plomax <- stan_model('dplomax.stan')
beta0_plomax <- vector()
beta1_plomax <- vector()

for(i in 1:r){
  set.seed(i)
  #400 samples
  fit_dplomax <- sampling(model_plomax, list(n = n, y=yy[,i], x=XX[,i]), iter = 200, chains = 4)
  
  ## Looic e Waic
  log_lik_1 <- extract_log_lik(fit_dplomax)
  measures[i,1] <- loo(log_lik_1)$looic
  measures[i,2] <- waic(log_lik_1)$waic
  
  #DIC, EAIC e EBIC
  postSamples <- Reduce(cbind, extract(fit_dplomax, pars=c("beta0", "beta1", "loglambda")))
  measures[i,3] <- dic(model.matrix(~XX[,i]), yy[,i], loglike_dplomax, postSamples)
  measures[i,4] <- eaic(model.matrix(~XX[,i]), yy[,i], loglike_dplomax, postSamples)
  measures[i,5] <- ebic(model.matrix(~XX[,i]), yy[,i], loglike_dplomax, postSamples)
  print(i)
}


write.csv(measures, "misspecification_dplomax25.csv")

##############################################################################
### Estimação Reversa de Potência Lomax ####################################################
###################################################################################
Finvplomax <- function(x, loglambda){
  lambda = exp(loglambda)
  p = 1 - F_dlomax(-x)^lambda
  return(p)
}

loglike_dpinvlomax <- function(par,y,x){
  eta <-x %*% par[-length(par)]
  p <- Finvplomax(eta, par[length(par)])
  logvero = sum(y*log(p) + (1-y)*log(1-p))
  return(logvero) 
}


model_rplomax <- stan_model('drplomax.stan')
beta0_rplomax <- vector()
beta1_rplomax <- vector()

for(i in 1:r){
  set.seed(i)
  #400 samples
  fit_rplomax <- sampling(model_rplomax, list(n = n, y=yy[,i], x=XX[,i]), iter = 200, chains = 4)
  
  ## Looic e Waic
  log_lik_1 <- extract_log_lik(fit_rplomax)
  measures[i,1] <- loo(log_lik_1)$looic
  measures[i,2] <- waic(log_lik_1)$waic
  
  #DIC, EAIC e EBIC
  postSamples <- Reduce(cbind, extract(fit_rplomax, pars=c("beta0", "beta1", "loglambda")))
  measures[i,3] <- dic(model.matrix(~XX[,i]), yy[,i] , loglike_dpinvlomax, postSamples)
  measures[i,4] <- eaic(model.matrix(~XX[,i]), yy[,i] , loglike_dpinvlomax, postSamples)
  measures[i,5] <- ebic(model.matrix(~XX[,i]), yy[,i] , loglike_dpinvlomax, postSamples)
  print(i)
}

write.csv(measures, "misspecification_dpinvlomax25.csv")

#########################################################################################
#### Analise dos resultados #######################################################
#################################################################################
### lambda = 2
mlog <- read.csv("misspecification_logistic4.csv")
mdlomax <- read.csv("misspecification_dlomax4.csv")
mdplomax <- read.csv("misspecification_dplomax4.csv")
mdpinvlomax <- read.csv("misspecification_dpinvlomax4.csv")

mean(mlog$DIC>mdlomax$DIC)
mean(mlog$EBIC>mdlomax$EBIC)
mean(mlog$EAIC>mdlomax$EAIC)
mean(mlog$WAIC > mdlomax$WAIC)
mean(mlog$LOO > mdlomax$LOO)

mean(mlog$DIC>mdplomax$DIC)
mean(mlog$EBIC>mdplomax$EBIC)
mean(mlog$EAIC>mdplomax$EAIC)
mean(mlog$WAIC > mdplomax$WAIC)
mean(mlog$LOO > mdplomax$LOO)

mean(mlog$DIC>mdpinvlomax$DIC)
mean(mlog$EBIC>mdpinvlomax$EBIC)
mean(mlog$EAIC>mdpinvlomax$EAIC)
mean(mlog$WAIC > mdpinvlomax$WAIC)
mean(mlog$LOO > mdpinvlomax$LOO)

mean(mlog$WAIC)
sd(mlog$WAIC)
mean(mdlomax$WAIC)
sd(mdlomax$WAIC)
mean(mdplomax$WAIC)
sd(mdplomax$WAIC)
mean(mdpinvlomax$WAIC)
sd(mdpinvlomax$WAIC)

mean(mlog$LOO)
sd(mlog$LOO)
mean(mdlomax$LOO)
sd(mdlomax$LOO)
mean(mdplomax$LOO)
sd(mdplomax$LOO)
mean(mdpinvlomax$LOO)
sd(mdpinvlomax$LOO)

