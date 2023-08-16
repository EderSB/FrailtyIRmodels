library(tidyverse)

#---------------------------------------------
# Estimation for the Frailty ARA_m model
#---------------------------------------------


#-------------------------------------------
# PLP functions
#-------------------------------------------

lambda.plp=function(t,beta,eta){
  resu=(beta/eta)*(t/eta)^(beta-1)
  return(resu)
}

Lambda.plp=function(t,beta,eta){
  resu=(t/eta)^beta
  return(resu)
}


#-------------------------------------------
# ARA_m functions
#-------------------------------------------

lambda.aram.plp=function(t,m,beta,eta,theta){
  time.aux=c(0,t)
  t.aux=c(rep(0,(m-1)),time.aux)
  v.age=NULL
  for(i in 1:length(t)){
    v.age[i]=t.aux[i+m]-(1-theta)*sum((theta^(0:(m-1)))*(rev(t.aux[i:(i+m-1)])))
  }
  value=lambda.plp(t=v.age, beta=beta, eta=eta)
  return(value)
}


Lambda.aram.plp<-function(t,m,beta,eta,theta){
  time.aux=c(0,t)
  t.aux=c(rep(0,(m-1)),time.aux)
  v.age=NULL
  v.age.aux=NULL
  resu.aux=NULL
  for(i in 1:length(t)){
    v.age[i]=t.aux[i+m]-(1-theta)*sum((theta^(0:(m-1)))*(rev(t.aux[i:(i+m-1)])))
    v.age.aux[i]=t.aux[i+m-1]-(1-theta)*sum((theta^(0:(m-1)))*(rev(t.aux[i:(i+m-1)])))
    resu.aux[i]=Lambda.plp(t=v.age[i],beta=beta,eta=eta)-
      Lambda.plp(t=v.age.aux[i],beta=beta,eta=eta)
  }
  resu=cumsum(resu.aux)
  return(resu)
}


#-------------------------------------------
# Frailty ARA_m functions
#-------------------------------------------

lambda.fr.aram.plp<-function(t,m,beta,eta,theta,alpha){
  value=lambda.aram.plp(t=t, m=m, beta=beta,eta=eta, theta=theta)/(1+alpha*Lambda.aram.plp(t=t,m=m, beta=beta,eta=eta, theta=theta))
  return(value)
}


Lambda.fr.aram.plp=function(t,m,beta,eta,theta,alpha){
  value = -log((1+alpha*Lambda.aram.plp(t=t,m=m,beta=beta,eta=eta,theta=theta))^(-1/alpha))
  return(value)
}



#------------------------------------------------------
# Log-likelihood function for the Frailty ARA_m model
#------------------------------------------------------

# Obs.: The "data" must be in data.frame, matrix or tibble format 
# containing at least three necessary information named by: 
# System (numbering of observed systems from 1 to the maximum 
# number of systems), Time (times until recorded failures) and 
# Event (1 for failure and 0 for truncation).
j=1
LogLik.aram <- function(par,data,m=1,par.fix=NULL){
  beta=exp(par[1])
  eta=exp(par[2])
  theta=1/(1+exp(-par[3]))
  alpha=exp(par[4])
  #
  k=max(data$System)
  l.k=NULL
  for(j in 1:k){
    data.k=data %>% filter(System==j)
    data.k.aux=data.k %>% filter(Event==1)
    time.fail=data.k.aux$Time
    time.obs=data.k$Time
    time.prev=c(0,time.obs[-length((time.obs))])
    
    lambda.aux=lambda.fr.aram.plp(t=time.fail,beta=beta,eta=eta,theta=theta,alpha=alpha,m=m)
    Lambda.aux=Lambda.fr.aram.plp(t=time.obs,beta=beta,eta=eta,theta=theta,alpha=alpha,m=m)- 
               Lambda.fr.aram.plp(t=time.prev,beta=beta,eta=eta,theta=theta,alpha=alpha,m=m)
   
    lambda.aux=ifelse(lambda.aux<0,1,lambda.aux)
    log.lambda.aux=log(lambda.aux)
    log.lambda.aux=ifelse(!is.finite(log.lambda.aux),0,log.lambda.aux)
    log.R.aux = -Lambda.aux
    #
    l.k[j] = -(sum(log.lambda.aux)+sum(log.R.aux)) 
  }
  l=sum(l.k)
  if(!is.finite(l)) l=1e50
  return(l %>% as.numeric())
  #
}


#------------------------------------------------------------------
# Log-likelihood optimization function for the Frailty ARA_m model
#------------------------------------------------------------------

#Obs.: We suggest the transformations c(log(beta),log(eta),-log(1/theta-1),log(alpha)) 
# for the initial.values

Analyse.aram <- function(data,m=NULL,initial.values,method='Nelder-Mead'){
  names=c('beta','eta','theta','alpha')
  npar=4
 
  mod=tryCatch(
    optim(
      fn=LogLik.aram, par=initial.values, hessian=TRUE,
      method=method, data=data, m=m),
    error=function(e) {e}
  )
  par=c(exp(mod$par[1]),exp(mod$par[2]),1/(1+exp(-mod$par[3])),exp(mod$par[4]))
  par=par %>% setNames(names)
  m2LL=-mod$value
  se0 = tryCatch(sqrt(diag(solve(mod$hessian, tol=1e-100))),
                 error=function(e){NA})
  #
  low_CI=exp(mod$par-1.96*se0)
  up_CI=exp(mod$par+1.96*se0)
  low_CI[3]=1/(1+exp(-mod$par[3]+1.96*se0[3]))
  up_CI[3]=1/(1+exp(-mod$par[3]-1.96*se0[3]))
  #
  AIC=-2*(m2LL)+2*length(par)
  BIC=-2*(m2LL)+length(par)*log(sum(data$Event))
  
  list(mod=mod, se0=se0, par.est=par, Conf.Int=round(data.frame(low_CI, up_CI),3), crit=m2LL, AIC=AIC, BIC=BIC)
  # ---
}


#------------------------------------------------------------------
# An example
#------------------------------------------------------------------

data = read.table("harvester_data.txt",h=TRUE)

result=Analyse.aram(data=data, m=3, initial.values=c(log(1.5),log(15),-log(1/0.4-1),log(0.2)),
                    method='Nelder-Mead')
result

