# clean the enviroment
rm(list=ls())

# load R packages
library("R2jags")
library("R2WinBUGS")
library(mvtnorm)

# set path of where the dm_gene is saved to source the function
setwd("/Users/qw/Desktop/Results_Final/cov8_vis2/cov8_lowpc_cor/")
# load data generation functions
source("dm_gene.R")

set.seed(666)


# parameter setup 
# we used 5 studies with random effects model

# number of studies
n.study<-5

# min and max sample size for each study
n_low<-100
n_up<-150


# total dimension of all fixed effects
# trt + cont + 8 baseline covariates
p.fix.all<-10

# dimension of baseline covariates
p.main<-8

# dimension of effect moderators
p.em<-8

# generate the design matrix
# generate design matrix
dm_input<-dm.generation(n.study=n.study,n_low=n_low,n_up=n_up,p.em=p.em,p.main=p.main)
data_dm<-dm_input$data_dm
data_dm_trt<-dm_input$data_dm_trt
data_dm_control<-dm_input$data_dm_control
n.size<-dm_input$n.size
n.total<-dm_input$n.total

# the number of true effect moderators
p.em.true<-2
# the number of false effect moderators
p.em.false<-6

# check whether the total number of effect moderators is correct
p.em==p.em.true+p.em.false

# set the true value of coefficients for baseline
beta<-c(2,3,1.8,2.7,2.3,1.5,1.7,2.2,1.3,2.6)

# set the true value of effect moderators
gamma.true<-rep(0.75,p.em.true)
gamma.false<-rep(0,p.em.false)

# put together coefficients for fixed effects
coef.study.fix<-c(beta,gamma.true,gamma.false)

# size of noise for all five studies
sigma_noise<-c(3.5,2.5,2.1,2.8,3.0)

# standard deviation of random effects for intercept, treatment, and effect moderators
tau_mu<-1.5
tau_trt<-1.5
tau_gamma<-c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5)

# dimension of random effects 
p.random<-10

# put together true values for all parameters
para.t<-c(coef.study.fix,tau_mu,tau_trt,tau_gamma,sigma_noise)

# generate outcomes for for each study
response<-list()
for (i in 1:n.study) {
  epsilon<-rnorm(n.size[i],0,sigma_noise[i])
  u_mu<-rnorm(1,0,tau_mu)
  u_trt<-rnorm(1,0,tau_trt)
  u_gamma<-rnorm(p.em,0,tau_gamma)
  coef.study<-c(coef.study.fix,u_mu,u_trt,u_gamma)
  response[[i]]<-as.matrix(data_dm[data_dm$study_id==i,-1])%*%coef.study+epsilon
}

# combine outcomes from five studies 
outcome<-do.call(rbind, response)

# put all datasets together
data_final<-data.frame(outcome,data_dm)

#############################################
#############################################
#############################################
##################  g=n for EM ##############
#############################################
#############################################
#############################################
model.gn = function(){
  
  # outcome model
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
  }
  
  # vague prior for main effects - coefficients of baseline covariates
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  # G=N prior for effect moderators
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,t(X[,(j+1)]*sqrt(phi[X[,1]]))%*%(X[,(j+1)]*sqrt(phi[X[,1]]))/n)
  }
  
  # distribution for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  # half-cauchy prior for the square root of variance component in the random effects model
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # vague gamma prior for the square root of precision
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-sqrt(1/phi[k])
  }
  
  #phi ~ dgamma(.005, .005)
  #sigma<-1/phi
}


  
# first column of data_final is outcome
data_gn = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study)

# implement Bayesian model with G=N   
gnout = jags(data=data_gn, model= model.gn,n.chain=2,n.burnin=10000,progress.bar="none",parameters.to.save=c("coef_fix","sigma","tau","re_u"),n.iter=20000)

# save posterior estimates
coef_fix_mean_gn<-gnout$BUGSoutput$summary[grep("coef_fix",rownames(gnout$BUGSoutput$summary)),"mean"]
coef_fix_median_gn<-gnout$BUGSoutput$summary[grep("coef_fix",rownames(gnout$BUGSoutput$summary)),"50%"]
coef_fix_lower_gn<-gnout$BUGSoutput$summary[grep("coef_fix",rownames(gnout$BUGSoutput$summary)),"2.5%"]
coef_fix_upper_gn<-gnout$BUGSoutput$summary[grep("coef_fix",rownames(gnout$BUGSoutput$summary)),"97.5%"]
  


#############################################
#############################################
#############################################
###############  g=sqrt(n) for EM ###########
#############################################
#############################################
#############################################

# model string
model.gsqrtn = function(){
  # outcome model
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
  }
  
  # vague prior for main effects - coefficients of baseline covariates
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  # g=sqrt(N) prior for effect moderators
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,t(X[,(j+1)]*sqrt(phi[X[,1]]))%*%(X[,(j+1)]*sqrt(phi[X[,1]]))/sqrt(n))
  }
  
  # for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  # half-cauchy prior for the square root of variance component in the random effects model
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # vague gamma prior for the square root of precision
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-sqrt(1/phi[k])
  }
}


# generate outcomes for each study
response<-list()
for (i in 1:n.study) {
    epsilon<-rnorm(n.size[i],0,sigma_noise[i])
    u_mu<-rnorm(1,0,tau_mu)
    u_trt<-rnorm(1,0,tau_trt)
    u_gamma<-rnorm(p.em,0,tau_gamma)
    coef.study<-c(coef.study.fix,u_mu,u_trt,u_gamma)
    response[[i]]<-as.matrix(data_dm[data_dm$study_id==i,-1])%*%coef.study+epsilon
}

# combine outcomes from all studies  
outcome<-do.call(rbind, response)

# put together all the data  
data_final<-data.frame(outcome,data_dm)

# prepare data list for MCMC 
data_gsqrtn = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study)
  
# implement Bayesian model with g=sqrt(n)
gsqrtnout = jags(data=data_gsqrtn, model= model.gsqrtn,n.chain=1,progress.bar="none",n.burnin=10000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"),n.iter=20000)

# save output
coef_fix_mean_gsqrtn<-gsqrtnout$BUGSoutput$summary[grep("coef_fix",rownames(gsqrtnout$BUGSoutput$summary)),"mean"]
coef_fix_median_gsqrtn<-gsqrtnout$BUGSoutput$summary[grep("coef_fix",rownames(gsqrtnout$BUGSoutput$summary)),"50%"]
coef_fix_lower_gsqrtn<-gsqrtnout$BUGSoutput$summary[grep("coef_fix",rownames(gsqrtnout$BUGSoutput$summary)),"2.5%"]
coef_fix_upper_gsqrtn<-gsqrtnout$BUGSoutput$summary[grep("coef_fix",rownames(gsqrtnout$BUGSoutput$summary)),"97.5%"]
  
 
########################################################
################  Zellner-Siow (ZS) prior for em #######
########################################################

model.zsem = function(){             
  # outcome model
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
  }
  
  # vague prior for main effects - coefficients of baseline covariates
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  # zs prior for effect moderators
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,t(X[,(j+1)]*sqrt(phi[X[,1]]))%*%(X[,(j+1)]*sqrt(phi[X[,1]]))*invgem[j-p.fix.all]/n)
  }
  
  # inver-gamma prior for parameter g
  for (k in 1:p.em) {
    invgem[k] ~ dgamma(.5, .5)
    gem[k]<-1/invgem[k]
    sf[k]<-gem[k]/(1+gem[k])
  }
  
  
  # distribution for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  # half-cauchy prior for the square root of variance component in the random effects model
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # vague gamma prior for the square root of precision
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-sqrt(1/phi[k])
  }
}


print(paste("zs only for em : ",j))

# generate outcomes for each study
response<-list()
for (i in 1:n.study) {
    epsilon<-rnorm(n.size[i],0,sigma_noise[i])
    u_mu<-rnorm(1,0,tau_mu)
    u_trt<-rnorm(1,0,tau_trt)
    u_gamma<-rnorm(p.em,0,tau_gamma)
    #u_gamma2<-rnorm(1,0,tau_gamma2)
    coef.study<-c(coef.study.fix,u_mu,u_trt,u_gamma)
    response[[i]]<-as.matrix(data_dm[data_dm$study_id==i,-1])%*%coef.study+epsilon
}

# combine outcomes from all studies  
outcome<-do.call(rbind, response)

# put together all the data  
data_final<-data.frame(outcome,data_dm)

# prepare data list for MCMC 
data_zsem = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study)

# implement Bayesian model with with ZS prior for parameter g
zsemout = jags(data=data_zsem, model= model.zsem,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"),n.iter=10000)
  
coef_fix.mean<-zsemout$BUGSoutput$summary[grep("coef_fix",rownames(zsemout$BUGSoutput$summary)),"mean"]
coef_fix.median<-zsemout$BUGSoutput$summary[grep("coef_fix",rownames(zsemout$BUGSoutput$summary)),"50%"]
coef_fix.lower<-zsemout$BUGSoutput$summary[grep("coef_fix",rownames(zsemout$BUGSoutput$summary)),"2.5%"]
coef_fix.upper<-zsemout$BUGSoutput$summary[grep("coef_fix",rownames(zsemout$BUGSoutput$summary)),"97.5%"]
  

#####################################################################
#####################################################################
#####################################################################
################  Calibrated ZS Prior for EM ########################
#####################################################################
#####################################################################
#####################################################################

# model string
model.czs = function(){

# outcome model
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
  }

# vague prior for main effects - coefficients of baseline covariates
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }

# uniform prior for sample size calibration
  # for fixed effects
  for (k in 1:n.study) {
    n_stu[k] ~ dunif(1,n.size[k])
  }
  n.vec<-rep(n_stu[1:n.study],n.size[1:n.study])
  
# calibrated ZS prior for EM  
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,t(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))%*%(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))*invgem[j-p.fix.all])
  }
  for (k in 1:p.em) {
    invgem[k] ~ dgamma(.5, .5)
    gem[k]<-1/invgem[k]
    sf[k]<-gem[k]/(1+gem[k])
  }

  # distribution for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
# half-cauchy prior for the square root of variance component in the random effects model
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
# vague gamma prior for the square root of precision
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-sqrt(1/phi[k])
  }
}

# prepare data list for MCMC 
data_czs = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study,n.size=n.size)

# implement Bayesian model with calibrated ZS prior
czsout = jags(data=data_czs, model= model.czs,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"),n.iter=10000)

# save results
coef_fix.mean<-czsout$BUGSoutput$summary[grep("coef_fix",rownames(czsout$BUGSoutput$summary)),"mean"]
coef_fix.median<-czsout$BUGSoutput$summary[grep("coef_fix",rownames(czsout$BUGSoutput$summary)),"50%"]
coef_fix.lower<-czsout$BUGSoutput$summary[grep("coef_fix",rownames(czsout$BUGSoutput$summary)),"2.5%"]
coef_fix.upper<-czsout$BUGSoutput$summary[grep("coef_fix",rownames(czsout$BUGSoutput$summary)),"97.5%"]
  

# implement Bayesian model with g=sqrt(n)

#####################################################################
#####################################################################
#####################################################################
########################  HG prior (a=3) ############################
#####################################################################
#####################################################################
#####################################################################

# model string
model_hg3 = function(){
  
  # outcome model
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
  }
  
  # vague prior for main effects - coefficients of baseline covariates
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  # HG(a=3) prior for EM
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,t(X[,(j+1)]*sqrt(phi[X[,1]]))%*%(X[,(j+1)]*sqrt(phi[X[,1]]))/gem[j-p.fix.all])
  }
  
  # prior for parameter g
  # replace a=3 with a=4 for HG (a=4)
  a = 3
  bw=a/2-1
  for (i in 1:p.em) {
    wem[i]~dbeta(1,bw)
    gem[i]<-wem[i]/(1-wem[i])
    sf[i]<-gem[i]/(1+gem[i])
  }
  
  # distribution for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  # half-cauchy prior for the square root of variance component in the random effects model
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # vague gamma prior for the square root of precision
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-sqrt(1/phi[k])
  }
}


# prepare data list for MCMC 
data_hg3 = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study)

# implement Bayesian model with HG prior
hg3out = jags(data=data_hg3, model=model_hg3,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"), 
                n.iter=10000)
  
coef_fix.mean.hg3<-hg3out$BUGSoutput$summary[grep("coef_fix",rownames(hg3out$BUGSoutput$summary)),"mean"]
coef_fix.median.hg3<-hg3out$BUGSoutput$summary[grep("coef_fix",rownames(hg3out$BUGSoutput$summary)),"50%"]
coef_fix.lower.hg3<-hg3out$BUGSoutput$summary[grep("coef_fix",rownames(hg3out$BUGSoutput$summary)),"2.5%"]
coef_fix.upper.hg3<-hg3out$BUGSoutput$summary[grep("coef_fix",rownames(hg3out$BUGSoutput$summary)),"97.5%"]
  

#####################################################################
#####################################################################
#####################################################################
########################  HGN prior (a=3) ###########################
#####################################################################
#####################################################################
#####################################################################

model_hgn3 = function(){
  
  # outcome model
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
  }
  
  # vague prior for main effects - coefficients of baseline covariates
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  # HGN(a=3) prior for EM
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,t(X[,(j+1)]*sqrt(phi[X[,1]]))%*%(X[,(j+1)]*sqrt(phi[X[,1]]))/gem[j-p.fix.all])
  }
  
  # change a=4 for HGN (a=4) prior
  a = 3
  bw=a/2-1
  
  for (i in 1:p.em) {
    wem[i]~dbeta(1,bw)
    gem[i]<-n*wem[i]/(1-wem[i])
    sf[i]<-gem[i]/(1+gem[i])
  }
  
  # half-cauchy prior for the square root of variance component in the random effects model
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  # half-cauchy prior for the square root of variance component in the random effects model
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # vague gamma prior for the square root of precision
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-sqrt(1/phi[k])
  }
}

# prepare data list for MCMC 
data_hgn3 = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study)

# implement Bayesian model with HGN prior (a=4)
hgn3out = jags(data=data_hgn3, model=model_hgn3,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"), 
                 n.iter=10000)
  
coef_fix.mean.hgn3<-hgn3out$BUGSoutput$summary[grep("coef_fix",rownames(hgn3out$BUGSoutput$summary)),"mean"]
coef_fix.median.hgn3<-hgn3out$BUGSoutput$summary[grep("coef_fix",rownames(hgn3out$BUGSoutput$summary)),"50%"]
coef_fix.lower.hgn3<-hgn3out$BUGSoutput$summary[grep("coef_fix",rownames(hgn3out$BUGSoutput$summary)),"2.5%"]
coef_fix.upper.hgn3<-hgn3out$BUGSoutput$summary[grep("coef_fix",rownames(hgn3out$BUGSoutput$summary)),"97.5%"]
  

#####################################################################
#####################################################################
#####################################################################
############################  CMG prior #############################
#####################################################################
#####################################################################
#####################################################################

model_hgshape1 = function(){
  
  # outcome model
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
  }
  
  # vague prior for main effects - coefficients of baseline covariates
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  
  for (i in 1:p.em) {
    b[i]~dunif(0,2)
    wem[i]~ dbeta(2,b[i])
    gem[i]<-wem[i]/(1-wem[i])
    sf[i]<-gem[i]/(1+gem[i])
  }
  
  # for fixed effects
  for (k in 1:n.study) {
    n_stu[k] ~ dunif(1,n.size[k])
  }
  n.vec<-rep(n_stu[1:n.study],n.size[1:n.study])
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,t(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))%*%(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))/gem[j-p.fix.all])
  }
  
  # distribution for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  # half-cauchy prior for the square root of variance component in the random effects model
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # vague gamma prior for the square root of precision
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-1/phi[k]
  }
}


  response<-list()
  for (i in 1:n.study) {
    epsilon<-rnorm(n.size[i],0,sigma_noise[i])
    u_mu<-rnorm(1,0,tau_mu)
    u_trt<-rnorm(1,0,tau_trt)
    u_gamma<-rnorm(p.em,0,tau_gamma)
    coef.study<-c(coef.study.fix,u_mu,u_trt,u_gamma)
    
    response[[i]]<-as.matrix(data_dm[data_dm$study_id==i,-1])%*%coef.study+epsilon
  }
  
  outcome<-do.call(rbind, response)
  data_final<-data.frame(outcome,data_dm)
  
  data_hgshape1 = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study,n.size=n.size)
  
  hgshape1out = jags(data=data_hgshape1, model=model_hgshape1,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"), 
                     n.iter=10000)
  
  coef_fix.mean.hgshape1<-hgshape1out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape1out$BUGSoutput$summary)),"mean"]
  coef_fix.median.hgshape1<-hgshape1out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape1out$BUGSoutput$summary)),"50%"]
  coef_fix.lower.hgshape1<-hgshape1out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape1out$BUGSoutput$summary)),"2.5%"]
  coef_fix.upper.hgshape1<-hgshape1out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape1out$BUGSoutput$summary)),"97.5%"]
  
  coef_fix_mean_hgshape1[j,]<-coef_fix.mean.hgshape1
  coef_fix_median_hgshape1[j,]<-coef_fix.median.hgshape1
  coef_fix_lower_hgshape1[j,]<-coef_fix.lower.hgshape1
  coef_fix_upper_hgshape1[j,]<-coef_fix.upper.hgshape1
  
  mse_hgshape1[j,]<-pmse_FE(model.name=hgshape1out,data_dm=data_dm,p.fix.all=p.fix.all,p.em=p.em,coef.study=coef.study,coef_fix_mean=coef_fix.mean.hgshape1,coef.study.fix=coef.study.fix,n.study=n.study)
  
  tau.mean.hgshape1<-hgshape1out$BUGSoutput$summary[grep("tau",rownames(hgshape1out$BUGSoutput$summary)),"mean"]
  tau.median.hgshape1<-hgshape1out$BUGSoutput$summary[grep("tau",rownames(hgshape1out$BUGSoutput$summary)),"50%"]
  tau.lower.hgshape1<-hgshape1out$BUGSoutput$summary[grep("tau",rownames(hgshape1out$BUGSoutput$summary)),"2.5%"]
  tau.upper.hgshape1<-hgshape1out$BUGSoutput$summary[grep("tau",rownames(hgshape1out$BUGSoutput$summary)),"97.5%"]
  
  tau_mean_hgshape1[j,]<-tau.mean.hgshape1
  tau_median_hgshape1[j,]<-tau.median.hgshape1
  tau_lower_hgshape1[j,]<-tau.lower.hgshape1
  tau_upper_hgshape1[j,]<-tau.upper.hgshape1
  
  
  
  sf.mean.hgshape1<-hgshape1out$BUGSoutput$summary[grep("sf",rownames(hgshape1out$BUGSoutput$summary)),"mean"]
  sf.median.hgshape1<-hgshape1out$BUGSoutput$summary[grep("sf",rownames(hgshape1out$BUGSoutput$summary)),"50%"]
  sf.lower.hgshape1<-hgshape1out$BUGSoutput$summary[grep("sf",rownames(hgshape1out$BUGSoutput$summary)),"2.5%"]
  sf.upper.hgshape1<-hgshape1out$BUGSoutput$summary[grep("sf",rownames(hgshape1out$BUGSoutput$summary)),"97.5%"]
  
  sf_mean_hgshape1[j,]<-sf.mean.hgshape1
  sf_median_hgshape1[j,]<-sf.median.hgshape1
  sf_lower_hgshape1[j,]<-sf.lower.hgshape1
  sf_upper_hgshape1[j,]<-sf.upper.hgshape1
  
}

Time1=Sys.time()
t_hgshape1<-Time1-Time0

para.mean.hgshape1<-cbind(coef_fix_mean_hgshape1,tau_mean_hgshape1,sf_mean_hgshape1)
para.median.hgshape1<-cbind(coef_fix_median_hgshape1,tau_median_hgshape1,sf_median_hgshape1)

para_lower_hgshape1<-cbind(coef_fix_lower_hgshape1,tau_lower_hgshape1)
para_upper_hgshape1<-cbind(coef_fix_upper_hgshape1,tau_upper_hgshape1)

para_hgshape1<-list(para.mean.hgshape1=para.mean.hgshape1,para_lower_hgshape1=para_lower_hgshape1,
                    para_upper_hgshape1=para_upper_hgshape1,mse_hgshape1=mse_hgshape1,
                    para.median.hgshape1=para.median.hgshape1,t_hgshape1=t_hgshape1)


########################################################
# adjusted g-type prior - shape uniform from 0 to 1 ####
########################################################

model_hgshape2 = function(){
  
  # sampling model for Y | beta, phi 
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
    #Y[i] ~ dnorm(X[i,1:6]%*%beta10+X[i,7:8]%*%beta11, phi)
  }
  
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  for (i in 1:p.em) {
    wem[i]~ dbeta(1,1)
    gem[i]<-wem[i]/(1-wem[i])
    sf[i]<-gem[i]/(1+gem[i])
  }
  
  # for fixed effects
  for (k in 1:n.study) {
    n_stu[k] ~ dunif(1,n.size[k])
  }
  n.vec<-rep(n_stu[1:n.study],n.size[1:n.study])
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,t(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))%*%(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))/gem[j-p.fix.all])
  }
  
  # for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # prior for phi to approximate p(phi) = 1/phi
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-1/phi[k]
  }
}

coef_fix_mean_hgshape2<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_median_hgshape2<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_lower_hgshape2<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_upper_hgshape2<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)

tau_mean_hgshape2<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_median_hgshape2<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_lower_hgshape2<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_upper_hgshape2<-matrix(NA,nrow = n.freq, ncol = p.random)

sf_mean_hgshape2<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_median_hgshape2<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_lower_hgshape2<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_upper_hgshape2<-matrix(NA,nrow = n.freq, ncol = p.em)

# sigma2_mean_hgshape2<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_median_hgshape2<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_lower_hgshape2<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_upper_hgshape2<-matrix(NA,nrow = n.freq, ncol = 1)

mse_hgshape2<-matrix(NA,nrow = n.freq, ncol = 4)


Time0=Sys.time()

for (j in 1:n.freq) {
  set.seed(j*6+9)
  print(paste("hgshape2: ",j))
  
  response<-list()
  for (i in 1:n.study) {
    epsilon<-rnorm(n.size[i],0,sigma_noise[i])
    u_mu<-rnorm(1,0,tau_mu)
    u_trt<-rnorm(1,0,tau_trt)
    u_gamma<-rnorm(p.em,0,tau_gamma)
    #u_gamma2<-rnorm(1,0,tau_gamma2)
    coef.study<-c(coef.study.fix,u_mu,u_trt,u_gamma)
    
    response[[i]]<-as.matrix(data_dm[data_dm$study_id==i,-1])%*%coef.study+epsilon
  }
  
  outcome<-do.call(rbind, response)
  data_final<-data.frame(outcome,data_dm)
  
  data_hgshape2 = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study,n.size=n.size)
  
  hgshape2out = jags(data=data_hgshape2, model=model_hgshape2,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"), 
                     n.iter=10000)
  
  coef_fix.mean.hgshape2<-hgshape2out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape2out$BUGSoutput$summary)),"mean"]
  coef_fix.median.hgshape2<-hgshape2out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape2out$BUGSoutput$summary)),"50%"]
  coef_fix.lower.hgshape2<-hgshape2out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape2out$BUGSoutput$summary)),"2.5%"]
  coef_fix.upper.hgshape2<-hgshape2out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape2out$BUGSoutput$summary)),"97.5%"]
  
  coef_fix_mean_hgshape2[j,]<-coef_fix.mean.hgshape2
  coef_fix_median_hgshape2[j,]<-coef_fix.median.hgshape2
  coef_fix_lower_hgshape2[j,]<-coef_fix.lower.hgshape2
  coef_fix_upper_hgshape2[j,]<-coef_fix.upper.hgshape2
  
  mse_hgshape2[j,]<-pmse_FE(model.name=hgshape2out,data_dm=data_dm,p.fix.all=p.fix.all,p.em=p.em,coef.study=coef.study,coef_fix_mean=coef_fix.mean.hgshape2,coef.study.fix=coef.study.fix,n.study=n.study)
  
  tau.mean.hgshape2<-hgshape2out$BUGSoutput$summary[grep("tau",rownames(hgshape2out$BUGSoutput$summary)),"mean"]
  tau.median.hgshape2<-hgshape2out$BUGSoutput$summary[grep("tau",rownames(hgshape2out$BUGSoutput$summary)),"50%"]
  tau.lower.hgshape2<-hgshape2out$BUGSoutput$summary[grep("tau",rownames(hgshape2out$BUGSoutput$summary)),"2.5%"]
  tau.upper.hgshape2<-hgshape2out$BUGSoutput$summary[grep("tau",rownames(hgshape2out$BUGSoutput$summary)),"97.5%"]
  
  tau_mean_hgshape2[j,]<-tau.mean.hgshape2
  tau_median_hgshape2[j,]<-tau.median.hgshape2
  tau_lower_hgshape2[j,]<-tau.lower.hgshape2
  tau_upper_hgshape2[j,]<-tau.upper.hgshape2
  
  sf.mean.hgshape2<-hgshape2out$BUGSoutput$summary[grep("sf",rownames(hgshape2out$BUGSoutput$summary)),"mean"]
  sf.median.hgshape2<-hgshape2out$BUGSoutput$summary[grep("sf",rownames(hgshape2out$BUGSoutput$summary)),"50%"]
  sf.lower.hgshape2<-hgshape2out$BUGSoutput$summary[grep("sf",rownames(hgshape2out$BUGSoutput$summary)),"2.5%"]
  sf.upper.hgshape2<-hgshape2out$BUGSoutput$summary[grep("sf",rownames(hgshape2out$BUGSoutput$summary)),"97.5%"]
  
  sf_mean_hgshape2[j,]<-sf.mean.hgshape2
  sf_median_hgshape2[j,]<-sf.median.hgshape2
  sf_lower_hgshape2[j,]<-sf.lower.hgshape2
  sf_upper_hgshape2[j,]<-sf.upper.hgshape2
  
  
}

Time1=Sys.time()
t_hgshape2<-Time1-Time0

para.mean.hgshape2<-cbind(coef_fix_mean_hgshape2,tau_mean_hgshape2,sf_mean_hgshape2)
para.median.hgshape2<-cbind(coef_fix_median_hgshape2,tau_median_hgshape2,sf_median_hgshape2)

para_lower_hgshape2<-cbind(coef_fix_lower_hgshape2,tau_lower_hgshape2,sf_lower_hgshape2)
para_upper_hgshape2<-cbind(coef_fix_upper_hgshape2,tau_upper_hgshape2,sf_upper_hgshape2)

para_hgshape2<-list(para.mean.hgshape2=para.mean.hgshape2,para_lower_hgshape2=para_lower_hgshape2,
                    para_upper_hgshape2=para_upper_hgshape2,mse_hgshape2=mse_hgshape2,
                    para.median.hgshape2=para.median.hgshape2,t_hgshape2=t_hgshape2)


############################################################
# adjusted g-type prior - shape with point mass around 0 ###
############################################################

model_hgshape3 = function(){
  
  # sampling model for Y | beta, phi 
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
    #Y[i] ~ dnorm(X[i,1:6]%*%beta10+X[i,7:8]%*%beta11, phi)
  }
  
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  for (i in 1:p.em) {
    b[i] ~ dunif(0,2)
    wem[i]~ dbeta(b[i],2)
    gem[i]<-wem[i]/(1-wem[i])
    sf[i]<-gem[i]/(1+gem[i])
  }
  
  # for fixed effects
  for (k in 1:n.study) {
    n_stu[k] ~ dunif(1,n.size[k])
  }
  n.vec<-rep(n_stu[1:n.study],n.size[1:n.study])
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,t(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))%*%(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))/gem[j-p.fix.all])
  }
  
  # for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # prior for phi to approximate p(phi) = 1/phi
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-1/phi[k]
  }
}

coef_fix_mean_hgshape3<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_median_hgshape3<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_lower_hgshape3<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_upper_hgshape3<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)

tau_mean_hgshape3<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_median_hgshape3<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_lower_hgshape3<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_upper_hgshape3<-matrix(NA,nrow = n.freq, ncol = p.random)

sf_mean_hgshape3<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_median_hgshape3<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_lower_hgshape3<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_upper_hgshape3<-matrix(NA,nrow = n.freq, ncol = p.em)

# sigma2_mean_hgshape3<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_median_hgshape3<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_lower_hgshape3<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_upper_hgshape3<-matrix(NA,nrow = n.freq, ncol = 1)

mse_hgshape3<-matrix(NA,nrow = n.freq, ncol = 4)


Time0=Sys.time()

for (j in 1:n.freq) {
  set.seed(j*6+9)
  print(paste("hgshape3: ",j))
  
  response<-list()
  for (i in 1:n.study) {
    epsilon<-rnorm(n.size[i],0,sigma_noise[i])
    u_mu<-rnorm(1,0,tau_mu)
    u_trt<-rnorm(1,0,tau_trt)
    u_gamma<-rnorm(p.em,0,tau_gamma)
    #u_gamma2<-rnorm(1,0,tau_gamma2)
    coef.study<-c(coef.study.fix,u_mu,u_trt,u_gamma)
    
    response[[i]]<-as.matrix(data_dm[data_dm$study_id==i,-1])%*%coef.study+epsilon
  }
  
  outcome<-do.call(rbind, response)
  data_final<-data.frame(outcome,data_dm)
  
  data_hgshape3 = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study,n.size=n.size)
  
  hgshape3out = jags(data=data_hgshape3, model=model_hgshape3,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"), 
                     n.iter=10000)
  
  coef_fix.mean.hgshape3<-hgshape3out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape3out$BUGSoutput$summary)),"mean"]
  coef_fix.median.hgshape3<-hgshape3out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape3out$BUGSoutput$summary)),"50%"]
  coef_fix.lower.hgshape3<-hgshape3out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape3out$BUGSoutput$summary)),"2.5%"]
  coef_fix.upper.hgshape3<-hgshape3out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape3out$BUGSoutput$summary)),"97.5%"]
  
  coef_fix_mean_hgshape3[j,]<-coef_fix.mean.hgshape3
  coef_fix_median_hgshape3[j,]<-coef_fix.median.hgshape3
  coef_fix_lower_hgshape3[j,]<-coef_fix.lower.hgshape3
  coef_fix_upper_hgshape3[j,]<-coef_fix.upper.hgshape3
  
  mse_hgshape3[j,]<-pmse_FE(model.name=hgshape3out,data_dm=data_dm,p.fix.all=p.fix.all,p.em=p.em,coef.study=coef.study,coef_fix_mean=coef_fix.mean.hgshape3,coef.study.fix=coef.study.fix,n.study=n.study)
  
  tau.mean.hgshape3<-hgshape3out$BUGSoutput$summary[grep("tau",rownames(hgshape3out$BUGSoutput$summary)),"mean"]
  tau.median.hgshape3<-hgshape3out$BUGSoutput$summary[grep("tau",rownames(hgshape3out$BUGSoutput$summary)),"50%"]
  tau.lower.hgshape3<-hgshape3out$BUGSoutput$summary[grep("tau",rownames(hgshape3out$BUGSoutput$summary)),"2.5%"]
  tau.upper.hgshape3<-hgshape3out$BUGSoutput$summary[grep("tau",rownames(hgshape3out$BUGSoutput$summary)),"97.5%"]
  
  tau_mean_hgshape3[j,]<-tau.mean.hgshape3
  tau_median_hgshape3[j,]<-tau.median.hgshape3
  tau_lower_hgshape3[j,]<-tau.lower.hgshape3
  tau_upper_hgshape3[j,]<-tau.upper.hgshape3
  
  
  sf.mean.hgshape3<-hgshape3out$BUGSoutput$summary[grep("sf",rownames(hgshape3out$BUGSoutput$summary)),"mean"]
  sf.median.hgshape3<-hgshape3out$BUGSoutput$summary[grep("sf",rownames(hgshape3out$BUGSoutput$summary)),"50%"]
  sf.lower.hgshape3<-hgshape3out$BUGSoutput$summary[grep("sf",rownames(hgshape3out$BUGSoutput$summary)),"2.5%"]
  sf.upper.hgshape3<-hgshape3out$BUGSoutput$summary[grep("sf",rownames(hgshape3out$BUGSoutput$summary)),"97.5%"]
  
  sf_mean_hgshape3[j,]<-sf.mean.hgshape3
  sf_median_hgshape3[j,]<-sf.median.hgshape3
  sf_lower_hgshape3[j,]<-sf.lower.hgshape3
  sf_upper_hgshape3[j,]<-sf.upper.hgshape3
  
  
}

Time1=Sys.time()
t_hgshape3<-Time1-Time0

para.mean.hgshape3<-cbind(coef_fix_mean_hgshape3,tau_mean_hgshape3,sf_mean_hgshape3)
para.median.hgshape3<-cbind(coef_fix_median_hgshape3,tau_median_hgshape3,sf_median_hgshape3)

para_lower_hgshape3<-cbind(coef_fix_lower_hgshape3,tau_lower_hgshape3,sf_lower_hgshape3)
para_upper_hgshape3<-cbind(coef_fix_upper_hgshape3,tau_upper_hgshape3,sf_upper_hgshape3)

para_hgshape3<-list(para.mean.hgshape3=para.mean.hgshape3,para_lower_hgshape3=para_lower_hgshape3,
                    para_upper_hgshape3=para_upper_hgshape3,mse_hgshape3=mse_hgshape3,
                    para.median.hgshape3=para.median.hgshape3,t_hgshape3=t_hgshape3)

############################################################
# adjusted g-type prior - uniform on SF & log size  ########
############################################################

model_hgshape4 = function(){
  
  # sampling model for Y | beta, phi 
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
    #Y[i] ~ dnorm(X[i,1:6]%*%beta10+X[i,7:8]%*%beta11, phi)
  }
  
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  for (i in 1:p.em) {
    #wem[i]~ dbeta(1,1)
    b[i] ~ dunif(0,2)
    wem[i]~ dbeta(2,b[i])
    gem[i]<-wem[i]/(1-wem[i])
    sf[i]<-gem[i]/(1+gem[i])
  }
  
  # for fixed effects
  for (k in 1:n.study) {
    n_stu_prep[k] ~ dunif(1,n.size[k])
    n_stu[k]<-log(n_stu_prep[k])
  }
  n.vec<-rep(n_stu[1:n.study],n.size[1:n.study])
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,t(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))%*%(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))/gem[j-p.fix.all])
  }
  
  # for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # prior for phi to approximate p(phi) = 1/phi
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-1/phi[k]
  }
}

coef_fix_mean_hgshape4<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_median_hgshape4<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_lower_hgshape4<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_upper_hgshape4<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)

tau_mean_hgshape4<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_median_hgshape4<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_lower_hgshape4<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_upper_hgshape4<-matrix(NA,nrow = n.freq, ncol = p.random)

sf_mean_hgshape4<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_median_hgshape4<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_lower_hgshape4<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_upper_hgshape4<-matrix(NA,nrow = n.freq, ncol = p.em)

# sigma2_mean_hgshape4<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_median_hgshape4<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_lower_hgshape4<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_upper_hgshape4<-matrix(NA,nrow = n.freq, ncol = 1)

mse_hgshape4<-matrix(NA,nrow = n.freq, ncol = 4)


Time0=Sys.time()

for (j in 1:n.freq) {
  set.seed(j*6+9)
  print(paste("hgshape4: ",j))
  
  response<-list()
  for (i in 1:n.study) {
    epsilon<-rnorm(n.size[i],0,sigma_noise[i])
    u_mu<-rnorm(1,0,tau_mu)
    u_trt<-rnorm(1,0,tau_trt)
    u_gamma<-rnorm(p.em,0,tau_gamma)
    #u_gamma2<-rnorm(1,0,tau_gamma2)
    coef.study<-c(coef.study.fix,u_mu,u_trt,u_gamma)
    
    response[[i]]<-as.matrix(data_dm[data_dm$study_id==i,-1])%*%coef.study+epsilon
  }
  
  outcome<-do.call(rbind, response)
  data_final<-data.frame(outcome,data_dm)
  
  data_hgshape4 = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study,n.size=n.size)
  
  hgshape4out = jags(data=data_hgshape4, model=model_hgshape4,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"), 
                     n.iter=10000)
  
  coef_fix.mean.hgshape4<-hgshape4out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape4out$BUGSoutput$summary)),"mean"]
  coef_fix.median.hgshape4<-hgshape4out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape4out$BUGSoutput$summary)),"50%"]
  coef_fix.lower.hgshape4<-hgshape4out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape4out$BUGSoutput$summary)),"2.5%"]
  coef_fix.upper.hgshape4<-hgshape4out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape4out$BUGSoutput$summary)),"97.5%"]
  
  coef_fix_mean_hgshape4[j,]<-coef_fix.mean.hgshape4
  coef_fix_median_hgshape4[j,]<-coef_fix.median.hgshape4
  coef_fix_lower_hgshape4[j,]<-coef_fix.lower.hgshape4
  coef_fix_upper_hgshape4[j,]<-coef_fix.upper.hgshape4
  
  mse_hgshape4[j,]<-pmse_FE(model.name=hgshape4out,data_dm=data_dm,p.fix.all=p.fix.all,p.em=p.em,coef.study=coef.study,coef_fix_mean=coef_fix.mean.hgshape4,coef.study.fix=coef.study.fix,n.study=n.study)
  
  tau.mean.hgshape4<-hgshape4out$BUGSoutput$summary[grep("tau",rownames(hgshape4out$BUGSoutput$summary)),"mean"]
  tau.median.hgshape4<-hgshape4out$BUGSoutput$summary[grep("tau",rownames(hgshape4out$BUGSoutput$summary)),"50%"]
  tau.lower.hgshape4<-hgshape4out$BUGSoutput$summary[grep("tau",rownames(hgshape4out$BUGSoutput$summary)),"2.5%"]
  tau.upper.hgshape4<-hgshape4out$BUGSoutput$summary[grep("tau",rownames(hgshape4out$BUGSoutput$summary)),"97.5%"]
  
  tau_mean_hgshape4[j,]<-tau.mean.hgshape4
  tau_median_hgshape4[j,]<-tau.median.hgshape4
  tau_lower_hgshape4[j,]<-tau.lower.hgshape4
  tau_upper_hgshape4[j,]<-tau.upper.hgshape4
  
  sf.mean.hgshape4<-hgshape4out$BUGSoutput$summary[grep("sf",rownames(hgshape4out$BUGSoutput$summary)),"mean"]
  sf.median.hgshape4<-hgshape4out$BUGSoutput$summary[grep("sf",rownames(hgshape4out$BUGSoutput$summary)),"50%"]
  sf.lower.hgshape4<-hgshape4out$BUGSoutput$summary[grep("sf",rownames(hgshape4out$BUGSoutput$summary)),"2.5%"]
  sf.upper.hgshape4<-hgshape4out$BUGSoutput$summary[grep("sf",rownames(hgshape4out$BUGSoutput$summary)),"97.5%"]
  
  sf_mean_hgshape4[j,]<-sf.mean.hgshape4
  sf_median_hgshape4[j,]<-sf.median.hgshape4
  sf_lower_hgshape4[j,]<-sf.lower.hgshape4
  sf_upper_hgshape4[j,]<-sf.upper.hgshape4
  
}

Time1=Sys.time()
t_hgshape4<-Time1-Time0

para.mean.hgshape4<-cbind(coef_fix_mean_hgshape4,tau_mean_hgshape4,sf_mean_hgshape4)
para.median.hgshape4<-cbind(coef_fix_median_hgshape4,tau_median_hgshape4,sf_median_hgshape4)

para_lower_hgshape4<-cbind(coef_fix_lower_hgshape4,tau_lower_hgshape4,sf_lower_hgshape4)
para_upper_hgshape4<-cbind(coef_fix_upper_hgshape4,tau_upper_hgshape4,sf_upper_hgshape4)

para_hgshape4<-list(para.mean.hgshape4=para.mean.hgshape4,para_lower_hgshape4=para_lower_hgshape4,
                    para_upper_hgshape4=para_upper_hgshape4,mse_hgshape4=mse_hgshape4,
                    para.median.hgshape4=para.median.hgshape4,t_hgshape4=t_hgshape4)


############################################################
# adjusted g-type prior - uniform on SF & log size  ########
############################################################

model_hgshape5 = function(){
  
  # sampling model for Y | beta, phi 
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
    #Y[i] ~ dnorm(X[i,1:6]%*%beta10+X[i,7:8]%*%beta11, phi)
  }
  
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  for (i in 1:p.em) {
    wem[i]~ dbeta(1,1)
    gem[i]<-wem[i]/(1-wem[i])
    sf[i]<-gem[i]/(1+gem[i])
  }
  
  # for fixed effects
  for (k in 1:n.study) {
    n_stu_prep[k] ~ dunif(1,n.size[k])
    n_stu[k]<-log(n_stu_prep[k])
  }
  n.vec<-rep(n_stu[1:n.study],n.size[1:n.study])
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,t(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))%*%(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))/gem[j-p.fix.all])
  }
  
  # for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # prior for phi to approximate p(phi) = 1/phi
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-1/phi[k]
  }
}

coef_fix_mean_hgshape5<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_median_hgshape5<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_lower_hgshape5<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_upper_hgshape5<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)

tau_mean_hgshape5<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_median_hgshape5<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_lower_hgshape5<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_upper_hgshape5<-matrix(NA,nrow = n.freq, ncol = p.random)

sf_mean_hgshape5<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_median_hgshape5<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_lower_hgshape5<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_upper_hgshape5<-matrix(NA,nrow = n.freq, ncol = p.em)

# sigma2_mean_hgshape5<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_median_hgshape5<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_lower_hgshape5<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_upper_hgshape5<-matrix(NA,nrow = n.freq, ncol = 1)

mse_hgshape5<-matrix(NA,nrow = n.freq, ncol = 4)

Time0=Sys.time()

for (j in 1:n.freq) {
  set.seed(j*6+9)
  print(paste("hgshape5: ",j))
  
  response<-list()
  for (i in 1:n.study) {
    epsilon<-rnorm(n.size[i],0,sigma_noise[i])
    u_mu<-rnorm(1,0,tau_mu)
    u_trt<-rnorm(1,0,tau_trt)
    u_gamma<-rnorm(p.em,0,tau_gamma)
    #u_gamma2<-rnorm(1,0,tau_gamma2)
    coef.study<-c(coef.study.fix,u_mu,u_trt,u_gamma)
    
    response[[i]]<-as.matrix(data_dm[data_dm$study_id==i,-1])%*%coef.study+epsilon
  }
  
  outcome<-do.call(rbind, response)
  data_final<-data.frame(outcome,data_dm)
  
  data_hgshape5 = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study,n.size=n.size)
  
  hgshape5out = jags(data=data_hgshape5, model=model_hgshape5,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"), 
                     n.iter=10000)
  
  coef_fix.mean.hgshape5<-hgshape5out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape5out$BUGSoutput$summary)),"mean"]
  coef_fix.median.hgshape5<-hgshape5out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape5out$BUGSoutput$summary)),"50%"]
  coef_fix.lower.hgshape5<-hgshape5out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape5out$BUGSoutput$summary)),"2.5%"]
  coef_fix.upper.hgshape5<-hgshape5out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape5out$BUGSoutput$summary)),"97.5%"]
  
  coef_fix_mean_hgshape5[j,]<-coef_fix.mean.hgshape5
  coef_fix_median_hgshape5[j,]<-coef_fix.median.hgshape5
  coef_fix_lower_hgshape5[j,]<-coef_fix.lower.hgshape5
  coef_fix_upper_hgshape5[j,]<-coef_fix.upper.hgshape5
  
  mse_hgshape5[j,]<-pmse_FE(model.name=hgshape5out,data_dm=data_dm,p.fix.all=p.fix.all,p.em=p.em,coef.study=coef.study,coef_fix_mean=coef_fix.mean.hgshape5,coef.study.fix=coef.study.fix,n.study=n.study)
  
  tau.mean.hgshape5<-hgshape5out$BUGSoutput$summary[grep("tau",rownames(hgshape5out$BUGSoutput$summary)),"mean"]
  tau.median.hgshape5<-hgshape5out$BUGSoutput$summary[grep("tau",rownames(hgshape5out$BUGSoutput$summary)),"50%"]
  tau.lower.hgshape5<-hgshape5out$BUGSoutput$summary[grep("tau",rownames(hgshape5out$BUGSoutput$summary)),"2.5%"]
  tau.upper.hgshape5<-hgshape5out$BUGSoutput$summary[grep("tau",rownames(hgshape5out$BUGSoutput$summary)),"97.5%"]
  
  tau_mean_hgshape5[j,]<-tau.mean.hgshape5
  tau_median_hgshape5[j,]<-tau.median.hgshape5
  tau_lower_hgshape5[j,]<-tau.lower.hgshape5
  tau_upper_hgshape5[j,]<-tau.upper.hgshape5
  
  
  sf.mean.hgshape5<-hgshape5out$BUGSoutput$summary[grep("sf",rownames(hgshape5out$BUGSoutput$summary)),"mean"]
  sf.median.hgshape5<-hgshape5out$BUGSoutput$summary[grep("sf",rownames(hgshape5out$BUGSoutput$summary)),"50%"]
  sf.lower.hgshape5<-hgshape5out$BUGSoutput$summary[grep("sf",rownames(hgshape5out$BUGSoutput$summary)),"2.5%"]
  sf.upper.hgshape5<-hgshape5out$BUGSoutput$summary[grep("sf",rownames(hgshape5out$BUGSoutput$summary)),"97.5%"]
  
  sf_mean_hgshape5[j,]<-sf.mean.hgshape5
  sf_median_hgshape5[j,]<-sf.median.hgshape5
  sf_lower_hgshape5[j,]<-sf.lower.hgshape5
  sf_upper_hgshape5[j,]<-sf.upper.hgshape5
  
  
}

Time1=Sys.time()
t_hgshape5<-Time1-Time0

para.mean.hgshape5<-cbind(coef_fix_mean_hgshape5,tau_mean_hgshape5,sf_mean_hgshape5)
para.median.hgshape5<-cbind(coef_fix_median_hgshape5,tau_median_hgshape5,sf_median_hgshape5)

para_lower_hgshape5<-cbind(coef_fix_lower_hgshape5,tau_lower_hgshape5,sf_lower_hgshape5)
para_upper_hgshape5<-cbind(coef_fix_upper_hgshape5,tau_upper_hgshape5,sf_upper_hgshape5)

para_hgshape5<-list(para.mean.hgshape5=para.mean.hgshape5,para_lower_hgshape5=para_lower_hgshape5,
                    para_upper_hgshape5=para_upper_hgshape5,mse_hgshape5=mse_hgshape5,
                    para.median.hgshape5=para.median.hgshape5,t_hgshape5=t_hgshape5)

# hgshape 6
model_hgshape6 = function(){
  
  # sampling model for Y | beta, phi 
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
    #Y[i] ~ dnorm(X[i,1:6]%*%beta10+X[i,7:8]%*%beta11, phi)
  }
  
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  for (i in 1:p.em) {
    b[i] ~ dunif(0,2)
    wem[i]~ dbeta(b[i],2)
    gem[i]<-wem[i]/(1-wem[i])
    sf[i]<-gem[i]/(1+gem[i])
  }
  
  # for fixed effects
  for (k in 1:n.study) {
    n_stu_prep[k] ~ dunif(1,n.size[k])
    n_stu[k]<-log(n_stu_prep[k])
  }
  
  n.vec<-rep(n_stu[1:n.study],n.size[1:n.study])
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,t(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))%*%(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))/gem[j-p.fix.all])
  }
  
  # for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # prior for phi to approximate p(phi) = 1/phi
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-1/phi[k]
  }
}

coef_fix_mean_hgshape6<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_median_hgshape6<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_lower_hgshape6<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_upper_hgshape6<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)

tau_mean_hgshape6<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_median_hgshape6<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_lower_hgshape6<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_upper_hgshape6<-matrix(NA,nrow = n.freq, ncol = p.random)


sf_mean_hgshape6<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_median_hgshape6<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_lower_hgshape6<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_upper_hgshape6<-matrix(NA,nrow = n.freq, ncol = p.em)

# sigma2_mean_hgshape6<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_median_hgshape6<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_lower_hgshape6<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_upper_hgshape6<-matrix(NA,nrow = n.freq, ncol = 1)

mse_hgshape6<-matrix(NA,nrow = n.freq, ncol = 4)

Time0=Sys.time()

for (j in 1:n.freq) {
  set.seed(j*6+9)
  print(paste("hgshape 6: ",j))
  
  response<-list()
  for (i in 1:n.study) {
    epsilon<-rnorm(n.size[i],0,sigma_noise[i])
    u_mu<-rnorm(1,0,tau_mu)
    u_trt<-rnorm(1,0,tau_trt)
    u_gamma<-rnorm(p.em,0,tau_gamma)
    #u_gamma2<-rnorm(1,0,tau_gamma2)
    coef.study<-c(coef.study.fix,u_mu,u_trt,u_gamma)
    
    response[[i]]<-as.matrix(data_dm[data_dm$study_id==i,-1])%*%coef.study+epsilon
  }
  
  outcome<-do.call(rbind, response)
  data_final<-data.frame(outcome,data_dm)
  
  data_hgshape6 = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study,n.size=n.size)
  
  hgshape6out = jags(data=data_hgshape6, model=model_hgshape6,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"), 
                     n.iter=10000)
  
  coef_fix.mean.hgshape6<-hgshape6out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape6out$BUGSoutput$summary)),"mean"]
  coef_fix.median.hgshape6<-hgshape6out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape6out$BUGSoutput$summary)),"50%"]
  coef_fix.lower.hgshape6<-hgshape6out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape6out$BUGSoutput$summary)),"2.5%"]
  coef_fix.upper.hgshape6<-hgshape6out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape6out$BUGSoutput$summary)),"97.5%"]
  
  coef_fix_mean_hgshape6[j,]<-coef_fix.mean.hgshape6
  coef_fix_median_hgshape6[j,]<-coef_fix.median.hgshape6
  coef_fix_lower_hgshape6[j,]<-coef_fix.lower.hgshape6
  coef_fix_upper_hgshape6[j,]<-coef_fix.upper.hgshape6
  
  mse_hgshape6[j,]<-pmse_FE(model.name=hgshape6out,data_dm=data_dm,p.fix.all=p.fix.all,p.em=p.em,coef.study=coef.study,coef_fix_mean=coef_fix.mean.hgshape6,coef.study.fix=coef.study.fix,n.study=n.study)
  
  tau.mean.hgshape6<-hgshape6out$BUGSoutput$summary[grep("tau",rownames(hgshape6out$BUGSoutput$summary)),"mean"]
  tau.median.hgshape6<-hgshape6out$BUGSoutput$summary[grep("tau",rownames(hgshape6out$BUGSoutput$summary)),"50%"]
  tau.lower.hgshape6<-hgshape6out$BUGSoutput$summary[grep("tau",rownames(hgshape6out$BUGSoutput$summary)),"2.5%"]
  tau.upper.hgshape6<-hgshape6out$BUGSoutput$summary[grep("tau",rownames(hgshape6out$BUGSoutput$summary)),"97.5%"]
  
  tau_mean_hgshape6[j,]<-tau.mean.hgshape6
  tau_median_hgshape6[j,]<-tau.median.hgshape6
  tau_lower_hgshape6[j,]<-tau.lower.hgshape6
  tau_upper_hgshape6[j,]<-tau.upper.hgshape6
  
  
  sf.mean.hgshape6<-hgshape6out$BUGSoutput$summary[grep("sf",rownames(hgshape6out$BUGSoutput$summary)),"mean"]
  sf.median.hgshape6<-hgshape6out$BUGSoutput$summary[grep("sf",rownames(hgshape6out$BUGSoutput$summary)),"50%"]
  sf.lower.hgshape6<-hgshape6out$BUGSoutput$summary[grep("sf",rownames(hgshape6out$BUGSoutput$summary)),"2.5%"]
  sf.upper.hgshape6<-hgshape6out$BUGSoutput$summary[grep("sf",rownames(hgshape6out$BUGSoutput$summary)),"97.5%"]
  
  sf_mean_hgshape6[j,]<-sf.mean.hgshape6
  sf_median_hgshape6[j,]<-sf.median.hgshape6
  sf_lower_hgshape6[j,]<-sf.lower.hgshape6
  sf_upper_hgshape6[j,]<-sf.upper.hgshape6
  
  
}

Time1=Sys.time()
t_hgshape6<-Time1-Time0

para.mean.hgshape6<-cbind(coef_fix_mean_hgshape6,tau_mean_hgshape6,sf_mean_hgshape6)
para.median.hgshape6<-cbind(coef_fix_median_hgshape6,tau_median_hgshape6,sf_median_hgshape6)

para_lower_hgshape6<-cbind(coef_fix_lower_hgshape6,tau_lower_hgshape6,sf_lower_hgshape6)
para_upper_hgshape6<-cbind(coef_fix_upper_hgshape6,tau_upper_hgshape6,sf_upper_hgshape6)

para_hgshape6<-list(para.mean.hgshape6=para.mean.hgshape6,para_lower_hgshape6=para_lower_hgshape6,
                    para_upper_hgshape6=para_upper_hgshape6,mse_hgshape6=mse_hgshape6,
                    para.median.hgshape6=para.median.hgshape6,t_hgshape6=t_hgshape6)


# hgshape 7
model_hgshape7 = function(){
  
  # sampling model for Y | beta, phi 
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
    #Y[i] ~ dnorm(X[i,1:6]%*%beta10+X[i,7:8]%*%beta11, phi)
  }
  
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  for (i in 1:p.em) {
    b[i] ~ dunif(0,2)
    wem[i]~ dbeta(2,b[i])
    gem[i]<-wem[i]/(1-wem[i])
    sf[i]<-gem[i]/(1+gem[i])
  }
  
  # for fixed effects
  for (k in 1:n.study) {
    alpha[k]~ dunif(0,1)
    n_stu[k] <- n.size[k]^alpha[k]
  }
  n.vec<-rep(n_stu[1:n.study],n.size[1:n.study])
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,t(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))%*%(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))/gem[j-p.fix.all])
  }
  
  # for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # prior for phi to approximate p(phi) = 1/phi
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-1/phi[k]
  }
}

coef_fix_mean_hgshape7<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_median_hgshape7<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_lower_hgshape7<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_upper_hgshape7<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)

tau_mean_hgshape7<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_median_hgshape7<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_lower_hgshape7<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_upper_hgshape7<-matrix(NA,nrow = n.freq, ncol = p.random)

sf_mean_hgshape7<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_median_hgshape7<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_lower_hgshape7<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_upper_hgshape7<-matrix(NA,nrow = n.freq, ncol = p.em)

# sigma2_mean_hgshape7<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_median_hgshape7<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_lower_hgshape7<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_upper_hgshape7<-matrix(NA,nrow = n.freq, ncol = 1)

mse_hgshape7<-matrix(NA,nrow = n.freq, ncol = 4)

Time0=Sys.time()

for (j in 1:n.freq) {
  set.seed(j*6+9)
  print(paste("hgshape 7: ",j))
  
  response<-list()
  for (i in 1:n.study) {
    epsilon<-rnorm(n.size[i],0,sigma_noise[i])
    u_mu<-rnorm(1,0,tau_mu)
    u_trt<-rnorm(1,0,tau_trt)
    u_gamma<-rnorm(p.em,0,tau_gamma)
    #u_gamma2<-rnorm(1,0,tau_gamma2)
    coef.study<-c(coef.study.fix,u_mu,u_trt,u_gamma)
    
    response[[i]]<-as.matrix(data_dm[data_dm$study_id==i,-1])%*%coef.study+epsilon
  }
  
  outcome<-do.call(rbind, response)
  data_final<-data.frame(outcome,data_dm)
  
  data_hgshape7 = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study,n.size=n.size)
  
  hgshape7out = jags(data=data_hgshape7, model=model_hgshape7,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"), 
                     n.iter=10000)
  
  coef_fix.mean.hgshape7<-hgshape7out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape7out$BUGSoutput$summary)),"mean"]
  coef_fix.median.hgshape7<-hgshape7out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape7out$BUGSoutput$summary)),"50%"]
  coef_fix.lower.hgshape7<-hgshape7out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape7out$BUGSoutput$summary)),"2.5%"]
  coef_fix.upper.hgshape7<-hgshape7out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape7out$BUGSoutput$summary)),"97.5%"]
  
  coef_fix_mean_hgshape7[j,]<-coef_fix.mean.hgshape7
  coef_fix_median_hgshape7[j,]<-coef_fix.median.hgshape7
  coef_fix_lower_hgshape7[j,]<-coef_fix.lower.hgshape7
  coef_fix_upper_hgshape7[j,]<-coef_fix.upper.hgshape7
  
  mse_hgshape7[j,]<-pmse_FE(model.name=hgshape7out,data_dm=data_dm,p.fix.all=p.fix.all,p.em=p.em,coef.study=coef.study,coef_fix_mean=coef_fix.mean.hgshape7,coef.study.fix=coef.study.fix,n.study=n.study)
  
  tau.mean.hgshape7<-hgshape7out$BUGSoutput$summary[grep("tau",rownames(hgshape7out$BUGSoutput$summary)),"mean"]
  tau.median.hgshape7<-hgshape7out$BUGSoutput$summary[grep("tau",rownames(hgshape7out$BUGSoutput$summary)),"50%"]
  tau.lower.hgshape7<-hgshape7out$BUGSoutput$summary[grep("tau",rownames(hgshape7out$BUGSoutput$summary)),"2.5%"]
  tau.upper.hgshape7<-hgshape7out$BUGSoutput$summary[grep("tau",rownames(hgshape7out$BUGSoutput$summary)),"97.5%"]
  
  tau_mean_hgshape7[j,]<-tau.mean.hgshape7
  tau_median_hgshape7[j,]<-tau.median.hgshape7
  tau_lower_hgshape7[j,]<-tau.lower.hgshape7
  tau_upper_hgshape7[j,]<-tau.upper.hgshape7
  
  sf.mean.hgshape7<-hgshape7out$BUGSoutput$summary[grep("sf",rownames(hgshape7out$BUGSoutput$summary)),"mean"]
  sf.median.hgshape7<-hgshape7out$BUGSoutput$summary[grep("sf",rownames(hgshape7out$BUGSoutput$summary)),"50%"]
  sf.lower.hgshape7<-hgshape7out$BUGSoutput$summary[grep("sf",rownames(hgshape7out$BUGSoutput$summary)),"2.5%"]
  sf.upper.hgshape7<-hgshape7out$BUGSoutput$summary[grep("sf",rownames(hgshape7out$BUGSoutput$summary)),"97.5%"]
  
  sf_mean_hgshape7[j,]<-sf.mean.hgshape7
  sf_median_hgshape7[j,]<-sf.median.hgshape7
  sf_lower_hgshape7[j,]<-sf.lower.hgshape7
  sf_upper_hgshape7[j,]<-sf.upper.hgshape7
  
  
}

Time1=Sys.time()
t_hgshape7<-Time1-Time0

para.mean.hgshape7<-cbind(coef_fix_mean_hgshape7,tau_mean_hgshape7,sf_mean_hgshape7)
para.median.hgshape7<-cbind(coef_fix_median_hgshape7,tau_median_hgshape7,sf_median_hgshape7)

para_lower_hgshape7<-cbind(coef_fix_lower_hgshape7,tau_lower_hgshape7)
para_upper_hgshape7<-cbind(coef_fix_upper_hgshape7,tau_upper_hgshape7)

para_hgshape7<-list(para.mean.hgshape7=para.mean.hgshape7,para_lower_hgshape7=para_lower_hgshape7,
                    para_upper_hgshape7=para_upper_hgshape7,mse_hgshape7=mse_hgshape7,
                    para.median.hgshape7=para.median.hgshape7,t_hgshape7=t_hgshape7)

# hgshape 8
model_hgshape8 = function(){
  
  # sampling model for Y | beta, phi 
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
    #Y[i] ~ dnorm(X[i,1:6]%*%beta10+X[i,7:8]%*%beta11, phi)
  }
  
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  for (i in 1:p.em) {
    #b[i] ~ dunif(0,2)
    wem[i]~ dbeta(1,1)
    gem[i]<-wem[i]/(1-wem[i])
    sf[i]<-gem[i]/(1+gem[i])
  }
  
  # for fixed effects
  for (k in 1:n.study) {
    alpha[k]~ dunif(0,1)
    n_stu[k] <- n.size[k]^alpha[k]
  }
  n.vec<-rep(n_stu[1:n.study],n.size[1:n.study])
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,t(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))%*%(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))/gem[j-p.fix.all])
  }
  
  # for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # prior for phi to approximate p(phi) = 1/phi
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-1/phi[k]
  }
}

coef_fix_mean_hgshape8<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_median_hgshape8<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_lower_hgshape8<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_upper_hgshape8<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)

tau_mean_hgshape8<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_median_hgshape8<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_lower_hgshape8<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_upper_hgshape8<-matrix(NA,nrow = n.freq, ncol = p.random)

sf_mean_hgshape8<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_median_hgshape8<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_lower_hgshape8<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_upper_hgshape8<-matrix(NA,nrow = n.freq, ncol = p.em)

# sigma2_mean_hgshape8<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_median_hgshape8<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_lower_hgshape8<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_upper_hgshape8<-matrix(NA,nrow = n.freq, ncol = 1)

mse_hgshape8<-matrix(NA,nrow = n.freq, ncol = 4)

Time0=Sys.time()

for (j in 1:n.freq) {
  set.seed(j*6+9)
  print(paste("hgshape 8: ",j))
  
  response<-list()
  for (i in 1:n.study) {
    epsilon<-rnorm(n.size[i],0,sigma_noise[i])
    u_mu<-rnorm(1,0,tau_mu)
    u_trt<-rnorm(1,0,tau_trt)
    u_gamma<-rnorm(p.em,0,tau_gamma)
    #u_gamma2<-rnorm(1,0,tau_gamma2)
    coef.study<-c(coef.study.fix,u_mu,u_trt,u_gamma)
    
    response[[i]]<-as.matrix(data_dm[data_dm$study_id==i,-1])%*%coef.study+epsilon
  }
  
  outcome<-do.call(rbind, response)
  data_final<-data.frame(outcome,data_dm)
  
  data_hgshape8 = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study,n.size=n.size)
  
  hgshape8out = jags(data=data_hgshape8, model=model_hgshape8,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"), 
                     n.iter=10000)
  
  coef_fix.mean.hgshape8<-hgshape8out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape8out$BUGSoutput$summary)),"mean"]
  coef_fix.median.hgshape8<-hgshape8out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape8out$BUGSoutput$summary)),"50%"]
  coef_fix.lower.hgshape8<-hgshape8out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape8out$BUGSoutput$summary)),"2.5%"]
  coef_fix.upper.hgshape8<-hgshape8out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape8out$BUGSoutput$summary)),"97.5%"]
  
  coef_fix_mean_hgshape8[j,]<-coef_fix.mean.hgshape8
  coef_fix_median_hgshape8[j,]<-coef_fix.median.hgshape8
  coef_fix_lower_hgshape8[j,]<-coef_fix.lower.hgshape8
  coef_fix_upper_hgshape8[j,]<-coef_fix.upper.hgshape8
  
  mse_hgshape8[j,]<-pmse_FE(model.name=hgshape8out,data_dm=data_dm,p.fix.all=p.fix.all,p.em=p.em,coef.study=coef.study,coef_fix_mean=coef_fix.mean.hgshape8,coef.study.fix=coef.study.fix,n.study=n.study)
  
  tau.mean.hgshape8<-hgshape8out$BUGSoutput$summary[grep("tau",rownames(hgshape8out$BUGSoutput$summary)),"mean"]
  tau.median.hgshape8<-hgshape8out$BUGSoutput$summary[grep("tau",rownames(hgshape8out$BUGSoutput$summary)),"50%"]
  tau.lower.hgshape8<-hgshape8out$BUGSoutput$summary[grep("tau",rownames(hgshape8out$BUGSoutput$summary)),"2.5%"]
  tau.upper.hgshape8<-hgshape8out$BUGSoutput$summary[grep("tau",rownames(hgshape8out$BUGSoutput$summary)),"97.5%"]
  
  tau_mean_hgshape8[j,]<-tau.mean.hgshape8
  tau_median_hgshape8[j,]<-tau.median.hgshape8
  tau_lower_hgshape8[j,]<-tau.lower.hgshape8
  tau_upper_hgshape8[j,]<-tau.upper.hgshape8
  
  sf.mean.hgshape8<-hgshape8out$BUGSoutput$summary[grep("sf",rownames(hgshape8out$BUGSoutput$summary)),"mean"]
  sf.median.hgshape8<-hgshape8out$BUGSoutput$summary[grep("sf",rownames(hgshape8out$BUGSoutput$summary)),"50%"]
  sf.lower.hgshape8<-hgshape8out$BUGSoutput$summary[grep("sf",rownames(hgshape8out$BUGSoutput$summary)),"2.5%"]
  sf.upper.hgshape8<-hgshape8out$BUGSoutput$summary[grep("sf",rownames(hgshape8out$BUGSoutput$summary)),"97.5%"]
  
  sf_mean_hgshape8[j,]<-sf.mean.hgshape8
  sf_median_hgshape8[j,]<-sf.median.hgshape8
  sf_lower_hgshape8[j,]<-sf.lower.hgshape8
  sf_upper_hgshape8[j,]<-sf.upper.hgshape8
  
  
}

Time1=Sys.time()
t_hgshape8<-Time1-Time0

para.mean.hgshape8<-cbind(coef_fix_mean_hgshape8,tau_mean_hgshape8,sf_mean_hgshape8)
para.median.hgshape8<-cbind(coef_fix_median_hgshape8,tau_median_hgshape8,sf_median_hgshape8)

para_lower_hgshape8<-cbind(coef_fix_lower_hgshape8,tau_lower_hgshape8)
para_upper_hgshape8<-cbind(coef_fix_upper_hgshape8,tau_upper_hgshape8)

para_hgshape8<-list(para.mean.hgshape8=para.mean.hgshape8,para_lower_hgshape8=para_lower_hgshape8,
                    para_upper_hgshape8=para_upper_hgshape8,mse_hgshape8=mse_hgshape8,
                    para.median.hgshape8=para.median.hgshape8,t_hgshape8=t_hgshape8)


# hgshape 9
model_hgshape9 = function(){
  
  # sampling model for Y | beta, phi 
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
    #Y[i] ~ dnorm(X[i,1:6]%*%beta10+X[i,7:8]%*%beta11, phi)
  }
  
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  for (i in 1:p.em) {
    b[i] ~ dunif(0,2)
    wem[i]~ dbeta(b[i],2)
    gem[i]<-wem[i]/(1-wem[i])
    sf[i]<-gem[i]/(1+gem[i])
  }
  
  # for fixed effects
  for (k in 1:n.study) {
    alpha[k]~ dunif(0,1)
    n_stu[k] <- n.size[k]^alpha[k]
  }
  n.vec<-rep(n_stu[1:n.study],n.size[1:n.study])
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,t(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))%*%(X[,(j+1)]*sqrt(phi[X[,1]])/sqrt(n.vec))/gem[j-p.fix.all])
  }
  
  # for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # prior for phi to approximate p(phi) = 1/phi
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-1/phi[k]
  }
}

coef_fix_mean_hgshape9<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_median_hgshape9<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_lower_hgshape9<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_upper_hgshape9<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)

tau_mean_hgshape9<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_median_hgshape9<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_lower_hgshape9<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_upper_hgshape9<-matrix(NA,nrow = n.freq, ncol = p.random)

sf_mean_hgshape9<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_median_hgshape9<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_lower_hgshape9<-matrix(NA,nrow = n.freq, ncol = p.em)
sf_upper_hgshape9<-matrix(NA,nrow = n.freq, ncol = p.em)

# sigma2_mean_hgshape9<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_median_hgshape9<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_lower_hgshape9<-matrix(NA,nrow = n.freq, ncol = 1)
# sigma2_upper_hgshape9<-matrix(NA,nrow = n.freq, ncol = 1)

mse_hgshape9<-matrix(NA,nrow = n.freq, ncol = 4)

Time0=Sys.time()

for (j in 1:n.freq) {
  set.seed(j*6+9)
  print(paste("hgshape 9: ",j))
  
  response<-list()
  for (i in 1:n.study) {
    epsilon<-rnorm(n.size[i],0,sigma_noise[i])
    u_mu<-rnorm(1,0,tau_mu)
    u_trt<-rnorm(1,0,tau_trt)
    u_gamma<-rnorm(p.em,0,tau_gamma)
    #u_gamma2<-rnorm(1,0,tau_gamma2)
    coef.study<-c(coef.study.fix,u_mu,u_trt,u_gamma)
    
    response[[i]]<-as.matrix(data_dm[data_dm$study_id==i,-1])%*%coef.study+epsilon
  }
  
  outcome<-do.call(rbind, response)
  data_final<-data.frame(outcome,data_dm)
  
  data_hgshape9 = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study,n.size=n.size)
  
  hgshape9out = jags(data=data_hgshape9, model=model_hgshape9,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"), 
                     n.iter=10000)
  
  coef_fix.mean.hgshape9<-hgshape9out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape9out$BUGSoutput$summary)),"mean"]
  coef_fix.median.hgshape9<-hgshape9out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape9out$BUGSoutput$summary)),"50%"]
  coef_fix.lower.hgshape9<-hgshape9out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape9out$BUGSoutput$summary)),"2.5%"]
  coef_fix.upper.hgshape9<-hgshape9out$BUGSoutput$summary[grep("coef_fix",rownames(hgshape9out$BUGSoutput$summary)),"97.5%"]
  
  coef_fix_mean_hgshape9[j,]<-coef_fix.mean.hgshape9
  coef_fix_median_hgshape9[j,]<-coef_fix.median.hgshape9
  coef_fix_lower_hgshape9[j,]<-coef_fix.lower.hgshape9
  coef_fix_upper_hgshape9[j,]<-coef_fix.upper.hgshape9
  
  mse_hgshape9[j,]<-pmse_FE(model.name=hgshape9out,data_dm=data_dm,p.fix.all=p.fix.all,p.em=p.em,coef.study=coef.study,coef_fix_mean=coef_fix.mean.hgshape9,coef.study.fix=coef.study.fix,n.study=n.study)
  
  tau.mean.hgshape9<-hgshape9out$BUGSoutput$summary[grep("tau",rownames(hgshape9out$BUGSoutput$summary)),"mean"]
  tau.median.hgshape9<-hgshape9out$BUGSoutput$summary[grep("tau",rownames(hgshape9out$BUGSoutput$summary)),"50%"]
  tau.lower.hgshape9<-hgshape9out$BUGSoutput$summary[grep("tau",rownames(hgshape9out$BUGSoutput$summary)),"2.5%"]
  tau.upper.hgshape9<-hgshape9out$BUGSoutput$summary[grep("tau",rownames(hgshape9out$BUGSoutput$summary)),"97.5%"]
  
  tau_mean_hgshape9[j,]<-tau.mean.hgshape9
  tau_median_hgshape9[j,]<-tau.median.hgshape9
  tau_lower_hgshape9[j,]<-tau.lower.hgshape9
  tau_upper_hgshape9[j,]<-tau.upper.hgshape9
  
  
  sf.mean.hgshape9<-hgshape9out$BUGSoutput$summary[grep("sf",rownames(hgshape9out$BUGSoutput$summary)),"mean"]
  sf.median.hgshape9<-hgshape9out$BUGSoutput$summary[grep("sf",rownames(hgshape9out$BUGSoutput$summary)),"50%"]
  sf.lower.hgshape9<-hgshape9out$BUGSoutput$summary[grep("sf",rownames(hgshape9out$BUGSoutput$summary)),"2.5%"]
  sf.upper.hgshape9<-hgshape9out$BUGSoutput$summary[grep("sf",rownames(hgshape9out$BUGSoutput$summary)),"97.5%"]
  
  sf_mean_hgshape9[j,]<-sf.mean.hgshape9
  sf_median_hgshape9[j,]<-sf.median.hgshape9
  sf_lower_hgshape9[j,]<-sf.lower.hgshape9
  sf_upper_hgshape9[j,]<-sf.upper.hgshape9
  
}

Time1=Sys.time()
t_hgshape9<-Time1-Time0

para.mean.hgshape9<-cbind(coef_fix_mean_hgshape9,tau_mean_hgshape9,sf_mean_hgshape9)
para.median.hgshape9<-cbind(coef_fix_median_hgshape9,tau_median_hgshape9,sf_median_hgshape9)

para_lower_hgshape9<-cbind(coef_fix_lower_hgshape9,tau_lower_hgshape9,sf_lower_hgshape9)
para_upper_hgshape9<-cbind(coef_fix_upper_hgshape9,tau_upper_hgshape9,sf_upper_hgshape9)

para_hgshape9<-list(para.mean.hgshape9=para.mean.hgshape9,para_lower_hgshape9=para_lower_hgshape9,
                    para_upper_hgshape9=para_upper_hgshape9,mse_hgshape9=mse_hgshape9,
                    para.median.hgshape9=para.median.hgshape9,t_hgshape9=t_hgshape9)


# regular approach 
model.const = function(){
  
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
  }
  
  for (j in 1:(p.fix.all+p.em)) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  # for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-sqrt(1/phi[k])
  }
}

coef_fix_mean_const<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_median_const<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_lower_const<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_upper_const<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)

tau_mean_const<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_median_const<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_lower_const<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_upper_const<-matrix(NA,nrow = n.freq, ncol = p.random)

mse_const<-matrix(NA,nrow = n.freq, ncol = 4)

Time0=Sys.time()

for (j in 1:n.freq) {
  set.seed(j*6+9)
  print(paste("iteration constant: ",j))
  
  ### random error 
  response<-list()
  for (i in 1:n.study) {
    epsilon<-rnorm(n.size[i],0,sigma_noise[i])
    u_mu<-rnorm(1,0,tau_mu)
    u_trt<-rnorm(1,0,tau_trt)
    u_gamma<-rnorm(p.em,0,tau_gamma)
    #u_gamma2<-rnorm(1,0,tau_gamma2)
    coef.study<-c(coef.study.fix,u_mu,u_trt,u_gamma)
    
    response[[i]]<-as.matrix(data_dm[data_dm$study_id==i,-1])%*%coef.study+epsilon
  }
  
  outcome<-do.call(rbind, response)
  
  data_final<-data.frame(outcome,data_dm)
  
  data_const = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study)
  
  constout = jags(data=data_const, model=model.const,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"), 
                  n.iter=10000)
  
  coef_fix.mean<-constout$BUGSoutput$summary[grep("coef_fix",rownames(constout$BUGSoutput$summary)),"mean"]
  coef_fix.median<-constout$BUGSoutput$summary[grep("coef_fix",rownames(constout$BUGSoutput$summary)),"50%"]
  coef_fix.lower<-constout$BUGSoutput$summary[grep("coef_fix",rownames(constout$BUGSoutput$summary)),"2.5%"]
  coef_fix.upper<-constout$BUGSoutput$summary[grep("coef_fix",rownames(constout$BUGSoutput$summary)),"97.5%"]
  
  coef_fix_mean_const[j,]<-coef_fix.mean
  coef_fix_median_const[j,]<-coef_fix.median
  coef_fix_lower_const[j,]<-coef_fix.lower
  coef_fix_upper_const[j,]<-coef_fix.upper
  
  mse_const[j,]<-pmse_FE(model.name=constout,data_dm=data_dm,p.fix.all=p.fix.all,p.em=p.em,coef.study=coef.study,coef_fix_mean=coef_fix.mean,coef.study.fix=coef.study.fix,n.study=n.study)
  
  tau.mean<-constout$BUGSoutput$summary[grep("tau",rownames(constout$BUGSoutput$summary)),"mean"]
  tau.median<-constout$BUGSoutput$summary[grep("tau",rownames(constout$BUGSoutput$summary)),"50%"]
  tau.lower<-constout$BUGSoutput$summary[grep("tau",rownames(constout$BUGSoutput$summary)),"2.5%"]
  tau.upper<-constout$BUGSoutput$summary[grep("tau",rownames(constout$BUGSoutput$summary)),"97.5%"]
  
  tau_mean_const[j,]<-tau.mean
  tau_median_const[j,]<-tau.median
  tau_lower_const[j,]<-tau.lower
  tau_upper_const[j,]<-tau.upper
  
}

Time1=Sys.time()
t_const<-Time1-Time0

para.mean.const<-cbind(coef_fix_mean_const,tau_mean_const)
para.median.const<-cbind(coef_fix_median_const,tau_median_const)
para_lower_const<-cbind(coef_fix_lower_const,tau_lower_const)
para_upper_const<-cbind(coef_fix_upper_const,tau_upper_const)

para_const<-list(para.mean.const=para.mean.const,para_lower_const=para_lower_const,
                 para_upper_const=para_upper_const,mse_const=mse_const,
                 para.median.const=para.median.const,t_const=t_const)

#############################################
# horseshoe prior #
#############################################

model_hs = function(){
  
  # sampling model for Y | beta, phi 
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
    #Y[i] ~ dnorm(X[i,1:6]%*%beta10+X[i,7:8]%*%beta11, phi)
  }
  
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  ## horseshoe prior distribution for penalization parameter
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    lambda[j-p.fix.all] ~ dt(0, 1, 1);T(0,)
    coef_fix[j] ~ dnorm(0,lambda[j-p.fix.all]^(-2)*eta^(-2))
  }
  eta ~ dt(0, 1, 1);T(0,)
  
  # for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # prior for phi to approximate p(phi) = 1/phi
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-sqrt(1/phi[k])
  }
}

coef_fix_mean_hs<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_median_hs<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_lower_hs<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_upper_hs<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)

tau_mean_hs<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_median_hs<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_lower_hs<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_upper_hs<-matrix(NA,nrow = n.freq, ncol = p.random)

mse_hs<-matrix(NA,nrow = n.freq, ncol = 4)


Time0=Sys.time()

for (j in 1:n.freq) {
  set.seed(j*6+9)
  print(paste("horseshoe iteration: ",j))
  
  response<-list()
  for (i in 1:n.study) {
    epsilon<-rnorm(n.size[i],0,sigma_noise[i])
    u_mu<-rnorm(1,0,tau_mu)
    u_trt<-rnorm(1,0,tau_trt)
    u_gamma<-rnorm(p.em,0,tau_gamma)
    #u_gamma2<-rnorm(1,0,tau_gamma2)
    coef.study<-c(coef.study.fix,u_mu,u_trt,u_gamma)
    
    response[[i]]<-as.matrix(data_dm[data_dm$study_id==i,-1])%*%coef.study+epsilon
  }
  
  outcome<-do.call(rbind, response)
  data_final<-data.frame(outcome,data_dm)
  
  data_hs = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study)
  
  hsout = jags(data=data_hs, model=model_hs,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"), 
               n.iter=10000)
  
  coef_fix.mean.hs<-hsout$BUGSoutput$summary[grep("coef_fix",rownames(hsout$BUGSoutput$summary)),"mean"]
  coef_fix.median.hs<-hsout$BUGSoutput$summary[grep("coef_fix",rownames(hsout$BUGSoutput$summary)),"50%"]
  coef_fix.lower.hs<-hsout$BUGSoutput$summary[grep("coef_fix",rownames(hsout$BUGSoutput$summary)),"2.5%"]
  coef_fix.upper.hs<-hsout$BUGSoutput$summary[grep("coef_fix",rownames(hsout$BUGSoutput$summary)),"97.5%"]
  
  coef_fix_mean_hs[j,]<-coef_fix.mean.hs
  coef_fix_median_hs[j,]<-coef_fix.median.hs
  coef_fix_lower_hs[j,]<-coef_fix.lower.hs
  coef_fix_upper_hs[j,]<-coef_fix.upper.hs
  
  mse_hs[j,]<-pmse_FE(model.name=hsout,data_dm=data_dm,p.fix.all=p.fix.all,p.em=p.em,coef.study=coef.study,coef_fix_mean=coef_fix.mean.hs,coef.study.fix=coef.study.fix,n.study=n.study)
  
  tau.mean.hs<-hsout$BUGSoutput$summary[grep("tau",rownames(hsout$BUGSoutput$summary)),"mean"]
  tau.median.hs<-hsout$BUGSoutput$summary[grep("tau",rownames(hsout$BUGSoutput$summary)),"50%"]
  tau.lower.hs<-hsout$BUGSoutput$summary[grep("tau",rownames(hsout$BUGSoutput$summary)),"2.5%"]
  tau.upper.hs<-hsout$BUGSoutput$summary[grep("tau",rownames(hsout$BUGSoutput$summary)),"97.5%"]
  
  tau_mean_hs[j,]<-tau.mean.hs
  tau_median_hs[j,]<-tau.median.hs
  tau_lower_hs[j,]<-tau.lower.hs
  tau_upper_hs[j,]<-tau.upper.hs
  
}

Time1=Sys.time()
t_hs<-Time1-Time0

para.mean.hs<-cbind(coef_fix_mean_hs,tau_mean_hs)
para.median.hs<-cbind(coef_fix_median_hs,tau_median_hs)

para_lower_hs<-cbind(coef_fix_lower_hs,tau_lower_hs)
para_upper_hs<-cbind(coef_fix_upper_hs,tau_upper_hs)

para_hs<-list(para.mean.hs=para.mean.hs,para_lower_hs=para_lower_hs,
              para_upper_hs=para_upper_hs,mse_hs=mse_hs,
              para.median.hs=para.median.hs,t_hs=t_hs)


#############################################
#################### ssvs1 ##################
#############################################

model_ssvs1 = function(){
  
  # sampling model for Y | beta, phi 
  for (i in 1:n) {
    Y[i]~ dnorm(X[i,2:(p.fix.all+p.em+1)]%*%coef_fix+X[i,(p.fix.all+p.em+2):(p.fix.all+p.em+p.random+2-1)]%*%re_u[1:p.random,X[i,1]], phi[X[i,1]])
    #Y[i] ~ dnorm(X[i,1:6]%*%beta10+X[i,7:8]%*%beta11, phi)
  }
  
  for (j in 1:p.fix.all) {
    coef_fix[j] ~ dnorm(0,0.001)
  }
  
  ## normal mixture distribution for penalization parameter
  
  for (j in (p.fix.all+1):(p.fix.all+p.em)) {
    index[j-p.fix.all] ~ dbern(0.5)
    latent[j-p.fix.all,1] ~ dnorm(0,1/(eta^2)) 
    latent[j-p.fix.all,2] ~ dnorm(0,1/(psi*eta^2))
    coef_fix[j]<-latent[j-p.fix.all,index[j-p.fix.all]+1]
  }
  psi<-100
  eta ~ dunif(0,5)
  
  # for random effects
  for (j in 1:p.random) {
    for (k in 1:n.study) {
      re_u[j,k] ~ dnorm(0,1/tau[j]^2)
    }
  }
  
  for (j in 1:p.random) {
    tau[j] ~ dt(0, 1, 1);T(0,)
  }
  
  # prior for phi to approximate p(phi) = 1/phi
  for (k in 1:n.study) {
    phi[k] ~ dgamma(.005, .005)
    sigma[k]<-sqrt(1/phi[k])
  }
}

coef_fix_mean_ssvs1<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_median_ssvs1<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_lower_ssvs1<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)
coef_fix_upper_ssvs1<-matrix(NA,nrow = n.freq, ncol = p.fix.all+p.em)

tau_mean_ssvs1<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_median_ssvs1<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_lower_ssvs1<-matrix(NA,nrow = n.freq, ncol = p.random)
tau_upper_ssvs1<-matrix(NA,nrow = n.freq, ncol = p.random)

mse_ssvs1<-matrix(NA,nrow = n.freq, ncol = 4)


Time0=Sys.time()

for (j in 1:n.freq) {
  set.seed(j*6+9)
  print(paste("ssvs1 iteration: ",j))
  
  response<-list()
  for (i in 1:n.study) {
    epsilon<-rnorm(n.size[i],0,sigma_noise[i])
    u_mu<-rnorm(1,0,tau_mu)
    u_trt<-rnorm(1,0,tau_trt)
    u_gamma<-rnorm(p.em,0,tau_gamma)
    #u_gamma2<-rnorm(1,0,tau_gamma2)
    coef.study<-c(coef.study.fix,u_mu,u_trt,u_gamma)
    
    response[[i]]<-as.matrix(data_dm[data_dm$study_id==i,-1])%*%coef.study+epsilon
  }
  
  outcome<-do.call(rbind, response)
  data_final<-data.frame(outcome,data_dm)
  
  data_ssvs1 = list(Y=data_final[,1], X=data_final[,-1],n=n.total,p.fix.all=p.fix.all,p.random=p.random,p.em=p.em,n.study=n.study)
  
  ssvs1out = jags(data=data_ssvs1, model=model_ssvs1,n.chain=1,n.burnin=5000,parameters.to.save=c("coef_fix","tau","sigma","re_u","sf"), 
                  n.iter=10000)
  
  coef_fix.mean.ssvs1<-ssvs1out$BUGSoutput$summary[grep("coef_fix",rownames(ssvs1out$BUGSoutput$summary)),"mean"]
  coef_fix.median.ssvs1<-ssvs1out$BUGSoutput$summary[grep("coef_fix",rownames(ssvs1out$BUGSoutput$summary)),"50%"]
  coef_fix.lower.ssvs1<-ssvs1out$BUGSoutput$summary[grep("coef_fix",rownames(ssvs1out$BUGSoutput$summary)),"2.5%"]
  coef_fix.upper.ssvs1<-ssvs1out$BUGSoutput$summary[grep("coef_fix",rownames(ssvs1out$BUGSoutput$summary)),"97.5%"]
  
  coef_fix_mean_ssvs1[j,]<-coef_fix.mean.ssvs1
  coef_fix_median_ssvs1[j,]<-coef_fix.median.ssvs1
  coef_fix_lower_ssvs1[j,]<-coef_fix.lower.ssvs1
  coef_fix_upper_ssvs1[j,]<-coef_fix.upper.ssvs1
  
  mse_ssvs1[j,]<-pmse_FE(model.name=ssvs1out,data_dm=data_dm,p.fix.all=p.fix.all,p.em=p.em,coef.study=coef.study,coef_fix_mean=coef_fix.mean.ssvs1,coef.study.fix=coef.study.fix,n.study=n.study)
  
  tau.mean.ssvs1<-ssvs1out$BUGSoutput$summary[grep("tau",rownames(ssvs1out$BUGSoutput$summary)),"mean"]
  tau.median.ssvs1<-ssvs1out$BUGSoutput$summary[grep("tau",rownames(ssvs1out$BUGSoutput$summary)),"50%"]
  tau.lower.ssvs1<-ssvs1out$BUGSoutput$summary[grep("tau",rownames(ssvs1out$BUGSoutput$summary)),"2.5%"]
  tau.upper.ssvs1<-ssvs1out$BUGSoutput$summary[grep("tau",rownames(ssvs1out$BUGSoutput$summary)),"97.5%"]
  
  tau_mean_ssvs1[j,]<-tau.mean.ssvs1
  tau_median_ssvs1[j,]<-tau.median.ssvs1
  tau_lower_ssvs1[j,]<-tau.lower.ssvs1
  tau_upper_ssvs1[j,]<-tau.upper.ssvs1
  
}

Time1=Sys.time()
t_ssvs1<-Time1-Time0

para.mean.ssvs1<-cbind(coef_fix_mean_ssvs1,tau_mean_ssvs1)
para.median.ssvs1<-cbind(coef_fix_median_ssvs1,tau_median_ssvs1)

para_lower_ssvs1<-cbind(coef_fix_lower_ssvs1,tau_lower_ssvs1)
para_upper_ssvs1<-cbind(coef_fix_upper_ssvs1,tau_upper_ssvs1)

para_ssvs1<-list(para.mean.ssvs1=para.mean.ssvs1,para_lower_ssvs1=para_lower_ssvs1,
                 para_upper_ssvs1=para_upper_ssvs1,mse_ssvs1=mse_ssvs1,
                 para.median.ssvs1=para.median.ssvs1,t_ssvs1=t_ssvs1)

#setwd("/Users/qw/Desktop/Moderation Project/sim/Multiple/July/NEW/Weak/5Stu/weak_XSD_Midtau")

cov8_c8<-list(para.t=para.t,para.true=para.true,para_gn=para_gn,para_gsqrtn=para_gsqrtn,
              para_zsem=para_zsem,para_zsemmod=para_zsemmod,
              para_hg3=para_hg3,para_hg4=para_hg4,
              para_hgn3=para_hgn3,para_hgn4=para_hgn4,
              para_hgshape1=para_hgshape1,para_hgshape2=para_hgshape2,para_hgshape3=para_hgshape3,
              para_hgshape4=para_hgshape4,para_hgshape5=para_hgshape5,para_hgshape6=para_hgshape6,
              para_hgshape7=para_hgshape7,para_hgshape8=para_hgshape8,para_hgshape9=para_hgshape9,
              para_const=para_const,para_hs=para_hs,para_ssvs1=para_ssvs1)

save(cov8_c8,file="cov8_c8.RData")




