rm(list=ls())

setwd("~/Documents/GitHub/lmoments_redux/monte_carlo")

mc.name = "gev"
paper_path = "~/"

quantile.function <- function(u, par){
  loc = par[1]
  scale = par[2]
  shape  = par[3]
  
  if(shape == 0)
    return(loc  - scale*log(-log(u))) else return(loc + scale*(1- (-log(u))^shape)/shape)
}

#Gradient of quantile function
grad.quantile.function <- function(u,par){
  loc = par[1]
  scale = par[2]
  shape  = par[3]
  
  if(shape==0)
    cbind(1, -log(-log(u)), (-(1/2))*scale*log(-log(u))^2)
    else    cbind(1, (1 - (-log(u))^shape)/shape, -((scale*(1 - (-log(u))^shape))/shape^2) - (scale*(-log(u))^shape*log(-log(u)))/shape)

}

hessian.quantile.function <- function(u,par){

    loc = par[1]
    scale = par[2]
    shape  = par[3]
    
    if(shape==0)
      cbind(c(0, 0, 0), c(0, 
                          0, (-(1/2))*log(-log(u))^2), c(0, (-(1/2))*log(-log(u))^2, (-(1/3))*
                                                           scale*log(-log(u))^3)) else cbind(c(0, 0, 0), c(0, 
                        0, -((1 - (-log(u))^shape)/shape^2) - ((-log(u))^shape*
                                                                 log(-log(u)))/
                          shape), c(0, -((1 - (-log(u))^shape)/shape^2) - ((-log(u))^shape*
                                                                             log(-log(u)))/shape, 
                                    (2*scale*(1 - (-log(u))^shape))/
                                      shape^3 + (2*scale*(-log(u))^shape*log(-log(u)))/
                                      shape^2 - (scale*(-log(u))^shape*log(-log(u))^2)/shape))
}

#Gradient of the quantile gradient function Q'(u|\theta)
grad.qdf <-function(u,par)
{
  loc = par[1]
  scale = par[2]
  shape  = par[3]
  
  if(shape==0)
    cbind(0, -(1/(u*log(u))), -((scale*log(-log(u)))/(u*log(u))))
  else cbind(0, (-log(u))^(-1 + shape)/u, (scale*(-log(u))^(-1 + shape)*log(-log(u)))/u)
}


density.function <- function(y, par){
  loc = par[1]
  scale = par[2]
  shape  = par[3]
  
  if(shape == 0)
    tY = exp(-(y-loc)/scale) else tY = (1- shape*(y-loc)/scale)^(1/shape)
  
  return(exp(-tY)*(tY^(-shape+1))/scale)
}

lmoment.analytic <- function(par,L)
{
  loc = par[1]
  scale = par[2]
  shape  = par[3]
  
  lvec = 0:(L-1)
  if(shape==0)
    const = -digamma(1) + log(1+lvec) else const = (1-((1+lvec)^(-shape))*gamma(1+shape))/shape
  return((1/(1+lvec))*(loc + scale*const))
}


lmoment.deriv.analytic <-function(par,L)
{
  loc = par[1]
  scale = par[2]
  shape  = par[3]
  
  lvec = 0:(L-1)
  
  if(shape==0)
    cbind(1/(1+lvec), (-digamma(1) + log(1+lvec))/(1+lvec), - scale*(6*digamma(1)^2 + pi^2 +12*(-digamma(1))*log(1+lvec) + 6*log(1+lvec)^2)/(12*(1+lvec)))
   else   cbind(1/(1+lvec), (1-((1+lvec)^(-shape))*gamma(1+shape))/(shape*(1+lvec)), -scale*(1-(1+lvec)^(-shape)*gamma(1+shape))/(shape^2*(1+lvec)) +scale*((1+lvec)^(-shape)*gamma(1+shape)*log(1+lvec)-(1+lvec)^(-shape)*gamma(1+shape)*digamma(1+shape))/(shape*(1+lvec)))
}

lmoment.hessian.analytic <-function(par, l){
  l=l-1
  r = par[2]
  k= par[3]
  if(k==0)
    cbind(c(0, 0, 0), c(0, 
                        0, -((6*(-digamma(1))^2 + pi^2 + 12*(-digamma(1))*log(1 + l) + 
                                6*log(1 + l)^2)/(12*(1 + l)))), 
          c(0, -((6*(-digamma(1))^2 + pi^2 + 12*(-digamma(1))*log(1 + l) + 
                    6*log(1 + l)^2)/(12*(1 + l))), 
            (r*(2*(-digamma(1))^3 + 
                  (-digamma(1))*pi^2 + (6*(-digamma(1))^2 + pi^2)*log(1 + l) + 
                  6*(-digamma(1))*log(1 + l)^2 + 2*log(1 + l)^3 - 
                  2*psi(2, 1)))/(6*(1 + l))))
    else cbind(c(0, 0, 0), c(0, 
                        0, -((1 - (1 + l)^(-k)*gamma(1 + k))/(
                          (k^2)*(1 + l))) + ((1 + l)^(-k)*gamma(1 + k)*log(1 + l) - ((1 + l)^-k)*gamma(1 + k)*digamma( 1 + k))/(
                            k*(1 + l))),c(0,-((1 - (1 + l)^(-k)*gamma(1 + k))/(
                              (k^2)*(1 + l))) + ((1 + l)^(-k)*gamma(1 + k)*log(1 + l) - ((1 + l)^-k)*gamma(1 + k)*digamma( 1 + k))/(
                                k*(1 + l)),(2*r*(1 - gamma(1 + k)/(1 + l)^k))/(k^3*(1 + l)) - (2*r*((gamma(1 + k)*log(1 + l))/(1 + l)^k - (gamma(1 + k)*digamma( 1 + k))/(1 + l)^k))/
            (k^2*(1 + l)) + (r*(-((gamma(1 + k)*log(1 + l)^2)/(1 + l)^k) + (2*gamma(1 + k)*log(1 + l)*digamma( 1 + k))/(1 + l)^k - 
                                  (gamma(1 + k)*digamma( 1 + k)^2)/(1 + l)^k - (gamma(1 + k)*psi(1, 1 + k))/(1 + l)^k))/(k*(1 + l))))  

}

true.par = c(0,1,-0.2)

Nreps = 5000
SampleSize = c(50,100,500)

max.L = 100
tau.seq = c(0.5,0.9,0.99,0.999)

source("aux/make_data.R")

source("aux/simulation.R")
source("aux/gen_results.R")

source("aux/gen_results_linear.R")

source("aux/simulation_selection.R")
source("aux/gen_results_selection.R")

source("aux/gen_correction.R")
source("aux/gen_results_coverage.R")

