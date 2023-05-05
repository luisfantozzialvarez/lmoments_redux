rm(list=ls())

setwd("~/Dropbox/Lmoments_redux/Codes")

mc.name = "gev"
paper_path = "~/Dropbox/Aplicativos/Overleaf/l_moments_redux"

julia_path = "/home/luisalvarez/julia-1.7.3/bin"

quantile.function <- function(u, par){
  loc = par[1]
  scale = par[2]
  shape  = par[3]
  
  if(shape == 0)
    return(loc  - scale*log(-log(u))) else return(loc + scale*(1- (-log(u))^shape)/shape)
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
    const = -digamma(1) + ((1+lvec)^(-shape))*log(1+lvec) else const = (1-((1+lvec)^(-shape))*gamma(1+shape))/shape
  return((1/(1+lvec))*(loc + scale*const))
}

true.par = c(0,1,-0.2)

Nreps = 500
SampleSize = c(50,100,500)

max.L = 100
tau.seq = c(0.5,0.9,0.99,0.999)

source("simulation.R")
source("gen_results.R")

source("simulation_selection.R")
source("gen_results_selection.R")

