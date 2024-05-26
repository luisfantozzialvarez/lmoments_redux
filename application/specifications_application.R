quantile.function.gpd <- function(u, par){
  loc = 0
  scale = par[1]
  shape  = par[2]
  
  if(shape==0)
    return(loc - scale*log(1-u)) else return(loc + scale*(1 - (1-u)^shape)/shape )
}

density.function.gpd <- function(y, par){
  loc = 0
  scale = par[1]
  shape  = par[2]
  
  if(shape==0)
    Tx = (y-loc)/scale  else Tx = (-1/shape)*log(1- shape*(y-loc)/scale)
  
  return((1/scale)*exp(-(1-shape)*Tx))
}

lmoment.analytic.gpd <- function(par,L)
{
  loc = 0
  scale = par[1]
  shape  = par[2]
  
  l_vec = 0:(L-1)
  
  if(shape==0)
    const = -(digamma(1) - digamma(2+l_vec)/gamma(2+l_vec)) else const = (1 - (gamma(1+shape)*gamma(2+l_vec))/gamma(2+shape+l_vec))/shape
  
  return((1/(1+l_vec))*(loc + scale*const))
}



quantile.function.gev <- function(u, par){
  loc = par[1]
  scale = par[2]
  shape  = par[3]
  
  if(shape == 0)
    return(loc  - scale*log(-log(u))) else return(loc + scale*(1- (-log(u))^shape)/shape)
}

density.function.gev <- function(y, par){
  loc = par[1]
  scale = par[2]
  shape  = par[3]
  
  if(shape == 0)
    tY = exp(-(y-loc)/scale) else tY = (1- shape*(y-loc)/scale)^(1/shape)
  
  return(exp(-tY)*(tY^(-shape+1))/scale)
}

lmoment.analytic.gev <- function(par,L)
{
  loc = par[1]
  scale = par[2]
  shape  = par[3]
  
  lvec = 0:(L-1)
  if(shape==0)
    const = -digamma(1) + ((1+lvec)^(-shape))*log(1+lvec) else const = (1-((1+lvec)^(-shape))*gamma(1+shape))/shape
  return((1/(1+lvec))*(loc + scale*const))
}





quantile.function.gev.mix <- function(u, par){
  prob = par[4]
  scale1 = par[3]
  scale2 = par[7]
  
  if(sign(scale1)*sign(scale2)==1)
    if(scale1<0)
      qts = ifelse(u<=prob,-quantile.function.gev(1-u/prob, par[1:3]),quantile.function.gev((u-prob)/(1-prob), par[5:7]))
    else  qts =   ifelse(u<=prob,quantile.function.gev(u/prob, par[1:3]),-quantile.function.gev(1 - (u-prob)/(1-prob), par[5:7]))

  if(sign(scale1)*sign(scale2)==-1)
    if(scale1<0)
      qts = ifelse(u<=prob,quantile.function.gev(u/prob, par[5:7]),quantile.function.gev((u-prob)/(1-prob), par[1:3]))
    else  qts = ifelse(u<=prob,quantile.function.gev(u/prob, par[1:3]),quantile.function.gev((u-prob)/(1-prob), par[5:7]))
 
  return(qts)

}


density.function.gev.mix <- function(y, par){
  prob = par[4]

  scale1 = par[3]
  scale2 = par[7]
  
  if(sign(scale1)*sign(scale2)==1)
    if(scale1<0)
    {
      p1 = density.function.gev(-y, par[1:3])
      
      p2 = density.function.gev(y, par[5:7])
      
    } else {
      p1 = density.function.gev(y, par[1:3])
      
      p2 = density.function.gev(-y, par[5:7])

    }
  
  if(sign(scale1)*sign(scale2)==-1)
    if(scale1<0)
    {
      p1 = density.function.gev(y, par[5:7])
      
      p2 = density.function.gev(y, par[1:3])
      
    } else {
      p1 = density.function.gev(y, par[1:3])
      
      p2 = density.function.gev(y, par[5:7])
    }

  dens = ifelse(is.nan(p1)&is.nan(p2), NA, prob*ifelse(is.nan(p1),0, p1) + (1-prob)*ifelse(is.nan(p2),0, p2))
  

  return(dens)

}


lmoment.analytic.gev.mix <- function(par,L)
{
  prob = par[4]
  scale1 = par[3]
  scale2 = par[7]
  
  if(sign(scale1)*sign(scale2)==1)
    if(scale1<0)
      return(lmoment.analytic.gev.mix.minus.plus(par,L)) else return(lmoment.analytic.gev.mix.plus.minus(par,L)) 
  
  if(sign(scale1)*sign(scale2)==-1)
    if(scale1<0)
      return(lmoment.analytic.gev.mix.same(c(par[5:7],par[4], par[1:3]),L)) else return(lmoment.analytic.gev.mix.same(par,L))
}



lmoment.analytic.gev.mix.same <- function(par,L)
{
  #rev=T
  
  loc1 = par[1]
  scale1 = par[2]
  shape1  = par[3]
  prob = par[4]
  loc2 = par[5]
  scale2 = par[6]
  shape2  = par[7]
  
  lvec = 0:(L-1)
  
  part1 = lmoment.analytic.gev(c(loc1,scale1,shape1), L)
  
  
  # part1f = rep(0,L)
  # 
  # for(j in lvec)
  #   part1f = part1f + choose(lvec, j)*part1[j+1]
  
  part1f = ((prob)^(lvec+1))*part1
  
  part2 = lmoment.analytic.gev(c(loc2,scale2,shape2), L)
  
  
  part2f = rep(0,L)
  for(j in lvec)
    part2f= part2f + choose(lvec, j)*part2[j+1]*((prob/(1-prob)))^(lvec-j)
  
  
  part2f = ((1-prob)^(lvec+1))*part2f
  
  return(part1f+part2f)
}


lmoment.analytic.gev.mix.minus.plus <- function(par,L)
{
  #rev=T

  loc1 = par[1]
  scale1 = par[2]
  shape1  = par[3]
  prob = par[4]
  loc2 = par[5]
  scale2 = par[6]
  shape2  = par[7]

  lvec = 0:(L-1)

  part1 = lmoment.analytic.gev(c(loc1,scale1,shape1), L)


  part1f = rep(0,L)

  for(j in lvec)
   part1f = part1f + choose(lvec, j)*((-1)^j)*part1[j+1]

  part1f = -((prob)^(lvec+1))*part1f

  part2 = lmoment.analytic.gev(c(loc2,scale2,shape2), L)


  part2f = rep(0,L)
  for(j in lvec)
    part2f= part2f + choose(lvec, j)*part2[j+1]*((prob/(1-prob)))^(lvec-j)


  part2f = ((1-prob)^(lvec+1))*part2f

  return(part1f+part2f)
}


lmoment.analytic.gev.mix.plus.minus <- function(par,L)
{
  #rev=T
  
  loc1 = par[1]
  scale1 = par[2]
  shape1  = par[3]
  prob = par[4]
  loc2 = par[5]
  scale2 = par[6]
  shape2  = par[7]
  
  lvec = 0:(L-1)
  
  part1 = lmoment.analytic.gev(c(loc1,scale1,shape1), L)
  
  
  # part1f = rep(0,L)
  # 
  # for(j in lvec)
  #   part1f = part1f + choose(lvec, j)*part1[j+1]
  
  
  part1f = ((prob)^(lvec+1))*part1
  
  part2 = lmoment.analytic.gev(c(loc2,scale2,shape2), L)
  
  
  part2f = rep(0,L)
  for(j in lvec)
    part2f= part2f + choose(lvec, j)*part2[j+1]*((-1)^j)*((1/(1-prob)))^(lvec-j)
  
  
  part2f = -((1-prob)^(lvec+1))*part2f
  
  return(part1f+part2f)
}



jacob.density.function.gev.mix <-  function(y, par){
    prob = par[4]
    
    scale1 = par[3]
    scale2 = par[7]
    
    if(sign(scale1)*sign(scale2)==1)
      if(scale1<0)
      {
        #p1 = ad_jacobian( function(x) density.function.gev(-y, x), par[1:3])
        
        #p2 =ad_jacobian( function(x) density.function.gev(y, x), par[5:7])
        
        p1 = matrix(0, nrow = length(y), ncol = length(par[1:3]))
        dens1 = density.function.gev(-y, par[1:3])
        case1 = !is.nan(dens1)
        dens1[!case1] = 0
        p1[case1,] =ad_jacobian( function(x) density.function.gev(-y[case1], x), par[1:3])
        
        
        p2 = matrix(0, nrow = length(y), ncol = length(par[5:7]))
        dens2 = density.function.gev(y, par[5:7])
        case2 = !is.nan(dens2)
        dens2[!case2] = 0
        p2[case2, ] = ad_jacobian( function(x) density.function.gev(y[case2], x), par[5:7])
        
        
        grd = cbind(prob*p1, dens1- dens2, (1-prob)*p2)
        
      } else {
        
        p1 = matrix(0, nrow = length(y), ncol = length(par[1:3]))
        dens1 = density.function.gev(y, par[1:3])
        case1 = !is.nan(dens1)
        dens1[!case1] = 0
        p1[case1,] =ad_jacobian( function(x) density.function.gev(y[case1], x), par[1:3])
        
        
        p2 = matrix(0, nrow = length(y), ncol = length(par[5:7]))
        dens2 = density.function.gev(-y, par[5:7])
        case2 = !is.nan(dens2)
        dens2[!case2] = 0
        p2[case2, ] = ad_jacobian( function(x) density.function.gev(-y[case2], x), par[5:7])
        
        
        #p1 = ad_jacobian( function(x) density.function.gev(y, x), par[1:3])
        
        #p2 =ad_jacobian( function(x) density.function.gev(-y, x), par[5:7])
        
        grd = cbind(prob*p1, dens1- dens2, (1-prob)*p2)

      }
    
    if(sign(scale1)*sign(scale2)==-1)
      if(scale1<0)
      {
        
        p1 = matrix(0, nrow = length(y), ncol = length(par[5:7]))
        dens1 = density.function.gev(y, par[5:7])
        case1 = !is.nan(dens1)
        dens1[!case1] = 0
        p1[case1,] =ad_jacobian( function(x) density.function.gev(y[case1], x), par[5:7])
        
        
        p2 = matrix(0, nrow = length(y), ncol = length(par[1:3]))
        dens2 = density.function.gev(y, par[1:3])
        case2 = !is.nan(dens2)
        dens2[!case2] = 0
        p2[case2, ] = ad_jacobian( function(x) density.function.gev(y[case2], x), par[1:3])
        
        
        

        grd = cbind((1-prob)*p2, dens1- dens2, prob*p1)
        

      } else {
        
        p1 = matrix(0, nrow = length(y), ncol = length(par[1:3]))
        dens1 = density.function.gev(y, par[1:3])
        case1 = !is.nan(dens1)
        dens1[!case1] = 0
        p1[case1,] =ad_jacobian( function(x) density.function.gev(y[case1], x), par[1:3])
        
        
        p2 = matrix(0, nrow = length(y), ncol = length(par[5:7]))
        dens2 = density.function.gev(y, par[5:7])
        case2 = !is.nan(dens2)
        dens2[!case2] = 0
        p2[case2, ] = ad_jacobian( function(x) density.function.gev(y[case2], x), par[5:7])
        
        
              
        grd = cbind(prob*p1, dens1- dens2, (1-prob)*p2)
        
        
      }
    

    
    return(grd)
    
  }



jacob.quantile.function.gev.mix <- function(u, par){
  prob = par[4]
  scale1 = par[3]
  scale2 = par[7]
  
  
  if(sign(scale1)*sign(scale2)==1)
    if(scale1<0)
      qts = ad_jacobian(function(x) c(-quantile.function.gev(1-u[u<=x[4]]/x[4], x[1:3]), quantile.function.gev((u[u>x[4]]-x[4])/(1-x[4]), x[5:7])), par)
  else  qts =   ad_jacobian(function(x) c(quantile.function.gev(u[u<=x[4]]/x[4],x[1:3]), -quantile.function.gev(1 - (u[u>x[4]]-x[4])/(1-x[4]),x[5:7])), par)
  
  if(sign(scale1)*sign(scale2)==-1)
    if(scale1<0)
      qts = ad_jacobian(function(x) c(quantile.function.gev(u[u<=x[4]]/x[4], x[5:7]), quantile.function.gev((u[u>x[4]]-x[4])/(1-x[4]),x[1:3])), par)
  else  qts = ad_jacobian(function(x) c(quantile.function.gev(u[u<=x[4]]/x[4],x[1:3]), quantile.function.gev((u[u>x[4]]-x[4])/(1-x[4]),x[5:7])), par)
  
  return(qts)
  
}

