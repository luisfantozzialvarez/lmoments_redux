rm(list=ls())

setwd("~/Dropbox/Lmoments_redux/Codes")

mc.name = "gpd"
paper_path = "~/Dropbox/Aplicativos/Overleaf/l_moments_redux"

julia_path = "/home/luisalvarez/julia-1.7.3/bin"

quantile.function <- function(u, par){
  loc = 0
  scale = par[1]
  shape  = par[2]

  if(shape==0)
    return(loc - scale*log(1-u)) else return(loc + scale*(1 - (1-u)^shape)/shape )
}

density.function <- function(y, par){
  loc = 0
  scale = par[1]
  shape  = par[2]
  
  if(shape==0)
    Tx = (y-loc)/scale  else Tx = (-1/shape)*log(1- shape*(y-loc)/scale)
  
  return((1/scale)*exp(-(1-shape)*Tx))
}

lmoment.analytic <- function(par,L)
{
  loc = 0
  scale = par[1]
  shape  = par[2]
  
  l_vec = 0:(L-1)
  
  if(shape==0)
    const = -(digamma(1) - digamma(2+l_vec)/gamma(2+l_vec)) else const = (1 - (gamma(1+shape)*gamma(2+l_vec))/gamma(2+shape+l_vec))/shape
  
  return((1/(1+l_vec))*(loc + scale*const))
}

true.par = c(1,-0.2)

Nreps = 500
SampleSize = c(50,100,500)

max.L=100
tau.seq = c(0.5,0.9,0.99,0.999)

source("simulation.R")
source("gen_results.R")

source("simulation_selection.R")
source("gen_results_selection.R")
