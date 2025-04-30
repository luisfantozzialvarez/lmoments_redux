quantile.function.gpd <- function(u, par){
  loc = 0
  scale = par[1]
  shape  = par[2]
  
  if(shape==0)
    return(loc - scale*log(1-u)) else return(loc + scale*(1 - (1-u)^shape)/shape )
}

#Gradient of quantile function
grad.quantile.function.gpd <- function(u,par){
  loc = 0
  scale = par[1]
  shape  = par[2]
  
  if(shape==0)
    cbind(1, -log(1 - u), (-(1/2))*scale*log(1 - u)^2)[,-1]
  else cbind(1, (1 - (1 - u)^shape)/shape, -((scale*(1 - (1 - u)^shape))/shape^2) - (scale*(1 - u)^shape*log(1 - u))/shape)[,-1]
  
}

hessian.quantile.function.gpd <- function(u,par){
  
  loc = 0
  scale = par[1]
  shape  = par[2]
  
  if(shape==0)
    cbind(c(0, 0, 0), c(0, 0, (-(1/2))*log(1 - u)^2), c(0, (-(1/2))*log(1 - u)^2, (-(1/3))*scale*log(1 - u)^3))[,-1][-1,]
  else cbind(c(0, 0, 0), c(0, 0, -((1 - (1 - u)^shape)/shape^2) - ((1 - u)^shape*log(1 - u))/shape), c(0, -((1 - (1 - u)^shape)/shape^2) - ((1 - u)^shape*log(1 - u))/shape, 
                                                                                                       (2*scale*(1 - (1 - u)^shape))/shape^3 + (2*scale*(1 - u)^shape*log(1 - u))/shape^2 - (scale*(1 - u)^shape*log(1 - u)^2)/shape))[,-1][-1,]
  
  
  
  
}

#Gradient of the quantile gradient function Q'(u|\theta)
grad.qdf.gpd <-function(u,par)
{
  loc = 0
  scale = par[1]
  shape  = par[2]
  
  if(shape==0)
    cbind(0, 1/(1 - u), (scale*log(1 - u))/(1 - u))[,-1] else     cbind(0, (1 - u)^(-1 + shape), scale*(1 - u)^(-1 + shape)*log(1 - u))[,-1]
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
    const = -(digamma(1) - digamma(2+l_vec)) else const = (1 - (gamma(1+shape)*gamma(2+l_vec))/gamma(2+shape+l_vec))/shape
  
  return((1/(1+l_vec))*(loc + scale*const))
}

lmoment.deriv.analytic.gpd <- function(par,L)
{
  m = 0
  r = par[1]
  k  = par[2]
  
  l = 0:(L-1)
  
  if(k==0)
    cbind(1/(1 + l), ((-digamma(1)) + digamma( 2 + l))/(1 + l), 
          -((r*gamma(
            1 + l)*(6*(-digamma(1))^2 + pi^2 + 
                      12*(-digamma(1))*digamma( 2 + l) + 6*digamma( 2 + l)^2 - 
                      6*psi(1, 2 + l)))/
              (12*gamma(2 + l))))[,-1]
  else 
    cbind(1/(1 + l), (1 - (gamma(1 + k)*gamma(2 + l))/
                        gamma(2 + k + l))/(k*(1 + l)), 
          -((r*(1 - (gamma(1 + k)*gamma(2 + l))/gamma(2 + k + l)))/(k^2*(1 + 
                                                                           l))) + 
            (r*(-((gamma(1 + k)*gamma(2 + l)*digamma( 1 + k))/
                    gamma(2 + k + l)) + (gamma(1 + k)*gamma(2 + l)*
                                           digamma( 2 + k + l))/
                  gamma(2 + k + l)))/(k*(1 + l)))[,-1]
  
}

lmoment.hessian.analytic.gpd <- function(par,l)
{
  m = 0
  r = par[1]
  k  = par[2]
  
  l = l-1
  
  if(k==0)
    mat = cbind(c(0, 0, 0), c(0, 
                              0, -((gamma(
                                1 + l)*(6*(-digamma(1))^2 + pi^2 + 
                                          12*(-digamma(1))*digamma( 2 + l) + 6*digamma( 2 + l)^2 - 
                                          6*psi(1, 2 + l)))/(12*gamma(2 + l)))), 
                c(0, -((gamma(
                  1 + l)*(6*(-digamma(1))^2 + pi^2 + 
                            12*(-digamma(1))*digamma( 2 + l) + 
                            6*digamma( 2 + l)^2 - 6*psi(1, 2 + l)))/
                    (12*gamma(2 + l))), (1/6)*
                    r*((1/gamma(2 + l))*(gamma(
                      1 + l)*(2*(-digamma(1))^3 + (-digamma(1))*pi^2 + 
                                6*(-digamma(1))*digamma( 2 + l)^2 + 
                                
                                digamma( 
                                  2 + l)*(6*(-digamma(1))^2 + pi^2 - 12*psi(1, 2 + l)) - 
                                6*(-digamma(1))*psi(1, 2 + l) - 2*psi(2, 1))) + 
                        (2*(digamma( 2 + l)^3 + 
                              3*digamma( 2 + l)*psi(1, 2 + l) + 
                              psi(2, 2 + l)))/(1 + l))))
  else mat = cbind(c(0, 0, 0), c(0, 
                                 0, -((1 - (gamma(1 + k)*gamma(2 + l))/
                                         gamma(2 + k + l))/(k^2*(1 + l))) + 
                                   (-((gamma(1 + k)*gamma(2 + l)*digamma( 1 + k))/
                                        gamma(2 + k + l)) + (gamma(1 + k)*gamma(2 + l)*
                                                               digamma( 2 + k + l))/
                                      
                                      gamma(2 + k + l))/(k*(1 + 
                                                              l))), c(0, -((1 - (gamma(1 + k)*gamma(2 + l))/
                                                                              gamma(2 + k + l))/(k^2*(1 + l))) + 
                                                                        (-((gamma(1 + k)*gamma(2 + l)*digamma( 1 + k))/
                                                                             gamma(2 + k + l)) + (gamma(1 + k)*gamma(2 + l)*
                                                                                                    digamma( 2 + k + l))/
                                                                           gamma(2 + k + l))/(k*(1 + l)), (2*
                                                                                                             r*(1 - (gamma(1 + k)*gamma(2 + l))/gamma(2 + k + l)))/(k^3*(1 + 
                                                                                                                                                                           l)) - 
                                                                        (2*
                                                                           r*(-((gamma(1 + k)*gamma(2 + l)*digamma( 1 + k))/
                                                                                  gamma(2 + k + l)) + (gamma(1 + k)*gamma(2 + l)*
                                                                                                         digamma( 2 + k + l))/
                                                                                gamma(2 + k + l)))/(k^2*(1 + 
                                                                                                           l)) + (1/(k*(1 + 
                                                                                                                          l)))*(r*(-((gamma(1 + k)*gamma(2 + l)*
                                                                                                                                        digamma( 1 + k)^2)/gamma(2 + k + l)) + 
                                                                                                                                     (2*gamma(1 + k)*gamma(2 + l)*digamma( 1 + k)*
                                                                                                                                        digamma( 2 + k + l))/gamma(2 + k + l) - 
                                                                                                                                     (gamma(1 + k)*gamma(2 + l)*digamma( 2 + k + l)^2)/
                                                                                                                                     gamma(2 + k + l) - (gamma(1 + k)*gamma(2 + l)*
                                                                                                                                                           psi(1, 1 + k))/
                                                                                                                                     
                                                                                                                                     gamma(2 + k + l) + (gamma(1 + k)*gamma(2 + l)*
                                                                                                                                                           psi(1, 2 + k + l))/gamma(2 + k + l)))))
  
  return(mat[-1,][,-1])
}


quantile.function.gev <- function(u, par){
  loc = par[1]
  scale = par[2]
  shape  = par[3]
  
  if(shape == 0)
    return(loc  - scale*log(-log(u))) else return(loc + scale*(1- (-log(u))^shape)/shape)
}

#Gradient of quantile function
grad.quantile.function.gev <- function(u,par){
  loc = par[1]
  scale = par[2]
  shape  = par[3]
  
  if(shape==0)
    cbind(1, -log(-log(u)), (-(1/2))*scale*log(-log(u))^2)
  else    cbind(1, (1 - (-log(u))^shape)/shape, -((scale*(1 - (-log(u))^shape))/shape^2) - (scale*(-log(u))^shape*log(-log(u)))/shape)
  
}

hessian.quantile.function.gev <- function(u,par){
  
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
grad.qdf.gev <-function(u,par)
{
  loc = par[1]
  scale = par[2]
  shape  = par[3]
  
  if(shape==0)
    cbind(0, -(1/(u*log(u))), -((scale*log(-log(u)))/(u*log(u))))
  else cbind(0, (-log(u))^(-1 + shape)/u, (scale*(-log(u))^(-1 + shape)*log(-log(u)))/u)
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
    const = -digamma(1) + log(1+lvec) else const = (1-((1+lvec)^(-shape))*gamma(1+shape))/shape
  return((1/(1+lvec))*(loc + scale*const))
}


lmoment.deriv.analytic.gev <-function(par,L)
{
  loc = par[1]
  scale = par[2]
  shape  = par[3]
  
  lvec = 0:(L-1)
  
  if(shape==0)
    cbind(1/(1+lvec), (-digamma(1) + log(1+lvec))/(1+lvec), - scale*(6*digamma(1)^2 + pi^2 +12*(-digamma(1))*log(1+lvec) + 6*log(1+lvec)^2)/(12*(1+lvec)))
  else   cbind(1/(1+lvec), (1-((1+lvec)^(-shape))*gamma(1+shape))/(shape*(1+lvec)), -scale*(1-(1+lvec)^(-shape)*gamma(1+shape))/(shape^2*(1+lvec)) +scale*((1+lvec)^(-shape)*gamma(1+shape)*log(1+lvec)-(1+lvec)^(-shape)*gamma(1+shape)*digamma(1+shape))/(shape*(1+lvec)))
}

lmoment.hessian.analytic.gev <-function(par, l){
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




quantile.function.gev.mix <- function(u, par){
  prob = par[4]
  shape1 = par[3]
  shape2 = par[7]
  
  if(sign(shape1)*sign(shape2)==1)
    if(shape1<0)
      qts = ifelse(u<=prob,-quantile.function.gev(1-u/prob, par[1:3]),quantile.function.gev((u-prob)/(1-prob), par[5:7]))
    else  qts =   ifelse(u<=prob,quantile.function.gev(u/prob, par[1:3]),-quantile.function.gev(1 - (u-prob)/(1-prob), par[5:7]))

  if(sign(shape1)*sign(shape2)==-1)
    if(shape1<0)
      qts = ifelse(u<=prob,-quantile.function.gev(1-u/prob, par[1:3]),-quantile.function.gev(1-(u-prob)/(1-prob), par[5:7]))
    else  qts = ifelse(u<=prob,quantile.function.gev(u/prob, par[1:3]),quantile.function.gev((u-prob)/(1-prob), par[5:7]))
 
  return(qts)

}



density.function.gev.mix <- function(y, par){
  prob = par[4]
  
  shape1 = par[3]
  shape2 = par[7]
  
  if(sign(shape1)*sign(shape2)==1)
    if(shape1<0)
    {
      p1 = density.function.gev(-y, par[1:3])
      
      p2 = density.function.gev(y, par[5:7])
      
    } else {
      p1 = density.function.gev(y, par[1:3])
      
      p2 = density.function.gev(-y, par[5:7])
      
    }
  
  if(sign(shape1)*sign(shape2)==-1)
    if(shape1<0)
    {
      p1 = density.function.gev(-y, par[1:3])
      
      p2 = density.function.gev(-y, par[5:7])
      
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
  shape1 = par[3]
  shape2 = par[7]
  
  if(sign(shape1)*sign(shape2)==1)
    if(shape1<0)
      return(lmoment.analytic.gev.mix.minus.plus(par,L)) else return(lmoment.analytic.gev.mix.plus.minus(par,L)) 
  
  if(sign(shape1)*sign(shape2)==-1)
    if(shape1<0)
      return(lmoment.analytic.gev.mix.minus.minus(par,L)) else return(lmoment.analytic.gev.mix.same(par,L))
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




lmoment.analytic.gev.mix.minus.minus <- function(par,L)
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
    part2f= part2f + choose(lvec, j)*part2[j+1]*((-1)^j)*((1/(1-prob)))^(lvec-j)
  
  
  part2f = -((1-prob)^(lvec+1))*part2f

  
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
  
  shape1 = par[3]
  shape2 = par[7]
  
  if(sign(shape1)*sign(shape2)==1)
    if(shape1<0)
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
  
  if(sign(shape1)*sign(shape2)==-1)
    if(shape1<0)
    {
      
      p1 = matrix(0, nrow = length(y), ncol = length(par[1:3]))
      dens1 = density.function.gev(-y, par[1:3])
      case1 = !is.nan(dens1)
      dens1[!case1] = 0
      p1[case1,] =ad_jacobian( function(x) density.function.gev(-y[case1], x), par[1:3])
      
      
      p2 = matrix(0, nrow = length(y), ncol = length(par[5:7]))
      dens2 = density.function.gev(-y, par[5:7])
      case2 = !is.nan(dens2)
      dens2[!case2] = 0
      p2[case2, ] = ad_jacobian( function(x) density.function.gev(-y[case2], x), par[5:7])
      
      
      
      
      grd = cbind(prob*p1, dens1- dens2, (1-prob)*p2)
      
      
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

lmoment.deriv.gev.mix <-function(par,L){
  prob = par[4]
  shape1 = par[3]
  shape2 = par[7]
  
  if(sign(shape1)*sign(shape2)==1)
    if(shape1<0)
      return(ad_jacobian(function(x) lmoment.analytic.gev.mix.minus.plus(x,L), par)) else return(ad_jacobian(function(x) lmoment.analytic.gev.mix.plus.minus(x,L), par)) 
  
  if(sign(shape1)*sign(shape2)==-1)
    if(shape1<0)
      return(ad_jacobian(function(x) lmoment.analytic.gev.mix.minus.minus(x,L), par)) else return(ad_jacobian(function(x) lmoment.analytic.gev.mix.same(x,L), par))
  
}

lmoment.deriv.gev.mix <-function(par,L){
  prob = par[4]
  shape1 = par[3]
  shape2 = par[7]
  
  if(sign(shape1)*sign(shape2)==1)
    if(shape1<0)
      return(ad_jacobian(function(x) lmoment.analytic.gev.mix.minus.plus(x,L), par)) else return(ad_jacobian(function(x) lmoment.analytic.gev.mix.plus.minus(x,L), par)) 
  
  if(sign(shape1)*sign(shape2)==-1)
    if(shape1<0)
      return(ad_jacobian(function(x) lmoment.analytic.gev.mix.minus.minus(x,L), par)) else return(ad_jacobian(function(x) lmoment.analytic.gev.mix.same(x,L), par))
  
}

jacob.quantile.function.gev.mix <- function(u, par){
  prob = par[4]
  shape1 = par[3]
  shape2 = par[7]
  
  
  if(sign(shape1)*sign(shape2)==1)
    if(shape1<0)
    {
      pt1 = rbind(-grad.quantile.function.gev(1-u[u<=par[4]]/par[4], par[1:3]), matrix(0, ncol=3, nrow = sum(u>par[4])))
      
      mm = c(-(u[u<=par[4]]/par[4]^2)/density.function.gev(quantile.function.gev(1-u[u<=par[4]]/par[4], par[1:3]),par[1:3]),
             (-1/(1-par[4])  +(u[u>par[4]]-par[4])/(1-par[4])^2)/density.function.gev(quantile.function.gev((u[u>par[4]]-par[4])/(1-par[4]), par[5:7]),par[5:7]))
      
      
      pt2 = rbind(matrix(0,ncol=3,nrow=sum(u<=par[4])), grad.quantile.function.gev((u[u>par[4]]-par[4])/(1-par[4]),par[5:7]))
      
    } else {
      pt1= rbind(grad.quantile.function.gev(u[u<=par[4]]/par[4], par[1:3]), matrix(0, ncol=3, nrow = sum(u>par[4])))
      mm = c(-((u[u<=par[4]])/par[4]^2)/density.function.gev(quantile.function.gev(u[u<=par[4]]/par[4], par[1:3]),par[1:3]),
             -((1-u[u>par[4]])/(1-par[4])^2 )/density.function.gev(quantile.function.gev((1- u[u>par[4]])/(1-par[4]),par[5:7]),par[5:7]))
      pt2 = rbind(matrix(0,ncol=3,nrow=sum(u<=par[4])), -grad.quantile.function.gev((1- u[u>par[4]])/(1-par[4]),par[5:7]))
      
    }
  
  if(sign(shape1)*sign(shape2)==-1)
    if(shape1<0)
    {
      pt1 = rbind(-grad.quantile.function.gev(1-u[u<=par[4]]/par[4], par[1:3]), matrix(0, ncol=3, nrow = sum(u>par[4])))
      mm = c(-((u[u<=par[4]])/par[4]^2)/density.function.gev(quantile.function.gev(1-u[u<=par[4]]/par[4], par[1:3]),par[1:3]),
             -((1-u[u>par[4]])/(1-par[4])^2 )/density.function.gev(quantile.function.gev((1- u[u>par[4]])/(1-par[4]),par[5:7]),par[5:7]))
      pt2 = rbind(matrix(0,ncol=3,nrow=sum(u<=par[4])), -grad.quantile.function.gev((1- u[u>par[4]])/(1-par[4]),par[5:7]))
      
      
    }else{
      pt1= rbind(grad.quantile.function.gev(u[u<=par[4]]/par[4], par[1:3]), matrix(0, ncol=3, nrow = sum(u>par[4])))
      mm = c(-((u[u<=par[4]])/par[4]^2)/density.function.gev(quantile.function.gev(u[u<=par[4]]/par[4], par[1:3]),par[1:3]),
             (-1/(1-par[4])  +(u[u>par[4]]-par[4])/(1-par[4])^2)/density.function.gev(quantile.function.gev((u[u>par[4]]-par[4])/(1-par[4]), par[5:7]),par[5:7]))
      
      pt2 = rbind(matrix(0,ncol=3,nrow=sum(u<=par[4])), grad.quantile.function.gev((u[u>par[4]]-par[4])/(1-par[4]),par[5:7]))
      
    }
  
  qts =cbind(pt1,mm,pt2)
  return(qts)
  
}
