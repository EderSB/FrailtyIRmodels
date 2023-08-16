library(tidyverse)

#---------------------------------------------
# Estimation for the Frailty ARI_m model
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
# ARI_m functions
#-------------------------------------------

lambda.arim.plp=function(t,m,beta,eta,theta){
  time.aux=c(0,t)
  lt=c(rep(0,(m-1)),lambda.plp(t=time.aux, beta=beta, eta=eta))
  value=NULL
  for(i in 1:length(t)){
    value[i]= lt[i+m]-(1-theta)*sum((theta^(0:(m-1)))*(rev(lt[i:(i+m-1)])))
  }
  return(value)
}


Lambda.arim.plp<-function(t,m,beta,eta,theta){
  time.aux=c(0,t)
  time.dif=diff(time.aux)
  time.aux2=time.aux[-length(time.aux)]
  lt=c(rep(0,(m-1)),lambda.plp(t=time.aux2, beta=beta, eta=eta))
  sum1=NULL
  for(i in 1:length(t)){
    sum1[i]=time.dif[i]*sum((theta^(0:(m-1)))*(rev(lt[i:(i+m-1)])))
  }
  resu=Lambda.plp(t=t,beta=beta,eta=eta)-(1-theta)*cumsum(sum1)
  return(resu)
}




#-------------------------------------------
# Frailty ARI_m functions
#-------------------------------------------

lambda.fr.arim.plp<-function(t,m,beta,eta,theta,alpha){
  value=lambda.arim.plp(t=t, m=m, beta=beta,eta=eta, theta=theta)/(1+alpha*Lambda.arim.plp(t=t,m=m, beta=beta,eta=eta, theta=theta))
  return(value)
}


Lambda.fr.arim.plp=function(t,m,beta,eta,theta,alpha){
  value = -log((1+alpha*Lambda.arim.plp(t=t,m=m,beta=beta,eta=eta,theta=theta))^(-1/alpha))
  return(value)
}


#------------------------------------------------------
# Log-likelihood function for the Frailty ARI_m model
#------------------------------------------------------

# Obs.: The "data" must be in data.frame, matrix or tibble format 
# containing at least three necessary information named by: 
# System (numbering of observed systems from 1 to the maximum 
# number of systems), Time (times until recorded failures) and 
# Event (1 for failure and 0 for truncation).

LogLik.arim <- function(par,data,m=1,par.fix=NULL){
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
    
    lambda.aux=lambda.fr.arim.plp(t=time.fail,beta=beta,eta=eta,theta=theta,alpha=alpha,m=m)
    Lambda.aux=Lambda.fr.arim.plp(t=time.obs,beta=beta,eta=eta,theta=theta,alpha=alpha,m=m)- 
               Lambda.fr.arim.plp(t=time.prev,beta=beta,eta=eta,theta=theta,alpha=alpha,m=m)
    
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
# Log-likelihood opimization function for the Frailty ARI_m model
#------------------------------------------------------------------

#Obs.: We suggest the transformations c(log(beta),log(eta),-log(1/theta-1),log(alpha)) 
# for the initial.values

Analyse.arim <- function(data,m=NULL,initial.values,method='Nelder-Mead'){
  names=c('beta','eta','theta','alpha')
  npar=4
  
  mod=tryCatch(
    optim(
      fn=LogLik.arim, par=initial.values, hessian=TRUE,
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
  
  list(mod=mod, se0=se0, par.plp=par, Conf.Int=round(data.frame(low_CI, up_CI),3), crit=m2LL, AIC=AIC, BIC=BIC)
  # ---
}


#------------------------------------------------------------------
# An example
#------------------------------------------------------------------

data = read.table("harvester_data.txt",h=TRUE)

result=Analyse.arim(data=data, m=3, initial.values=c(log(1.5),log(15),-log(1/0.4-1),log(0.2)),
                    method='Nelder-Mead')
result

