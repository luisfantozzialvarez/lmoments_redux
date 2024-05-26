library(np)
library(pracma)


lmoment.est <- function(y, par, L, orthogonal = F, lmoment.analytic = NULL, quantile.func=NULL,
                        weight.matrix = "par", density.function = NULL,  
                        lmoment.est = "caglad", grid.length = length(y), par.first.step = NULL, vcov= F,...){
 
  if(orthogonal)
    coefs =  sapply(0:(L-1),function(x){ c((-1)^(x-0:x)*sapply(0:x, function(y){exp(log(choose(x,y))+log(choose(x+y,y)) - log(x+1)*(x+1))}),
                                           rep(0, L-1-x))})
    
  y = y[order(y)]
  N = length(y)
  
  if(lmoment.est=="caglad")
    values = t(sapply(0:(L-1), function(x){
      vals = (((0:N)/N)^(x+1))/(x+1)
      sum(diff(vals)*y)
    })) else if(lmoment.est == "unbiased") values = y%*%sapply(0:(L-1), function(k){sapply(1:N, function(i){exp(log(choose(i-1,k))-log(choose(N,k+1)))})/(k+1)})
  
  
  if(orthogonal)
    l.est = values%*%coefs else l.est = values
  
  if(is.null(lmoment.analytic))
  {
   
    midpointer = (0:(grid.length-1) + 1:(grid.length))/(2*grid.length)
    
    powers = sapply(0:(L-1), function(x){
      midpointer^x})

    lmoment.func <- function(theta){
      vvv = colMeans(sapply(midpointer, function(u){quantile.func(u,theta)})*powers)
      if(orthogonal)
        return(vvv%*%coefs) else return(vvv)
    }
  } else {
    lmoment.func <- function(theta){
      lmoment.analytic(theta, L)
    }
  }
  
  objective <- function(theta){(lmoment.func(theta)-l.est)%*%(t(lmoment.func(theta)  -l.est))}
  
  if(is.null(par.first.step))
    first.step = optim(par, objective, method = "BFGS", ...) else first.step = list("par" = par.first.step)
  
  if(weight.matrix == "first")
    return(list("fs"=first.step)) else {
      if(weight.matrix == "semi")
      {
        grid.length.weight = N
        midder = (0:(grid.length.weight-1) + 1:(grid.length.weight))/(2*grid.length.weight)
        vwd = 1/sapply(y, function(x){density.function(x, first.step$par)})
        
      } else if(weight.matrix == "par")
      {
        grid.length.weight = grid.length
        midder = (0:(grid.length.weight-1) + 1:(grid.length.weight))/(2*grid.length.weight)
        vwd = sapply(midder, function(u){1/density.function(quantile.func(u,first.step$par), first.step$par)})
        
      } else if(weight.matrix == "np")
      {
        grid.length.weight = N
        midder = (0:(grid.length.weight-1) + 1:(grid.length.weight))/(2*grid.length.weight)
        
        bw = npudensbw(y)
        vwd = 1/fitted(npudens(bw))
      }
      
      
      mat_middle = outer(midder, rep(1,grid.length.weight))
      
      mat_middle = pmin(mat_middle, t(mat_middle)) - mat_middle*t(mat_middle)
      
      
      mat_middle = mat_middle*(vwd%*%t(vwd))
      
      power_mat = sapply(0:(L-1), function(l){midder^l})
      
      weight = t(power_mat)%*%mat_middle%*%power_mat/grid.length.weight^2
      
      if(orthogonal)
        weight = t(coefs)%*%weight%*%coefs
      
      weight = pinv(weight)
    }
    
  objective <- function(theta){(lmoment.func(theta)-l.est)%*%weight%*%(t(lmoment.func(theta)  -l.est))}
  
  second.step = optim(par, objective, method = "BFGS", ...)
  
  if(vcov)
  {
    hss = - ad_jacobian(lmoment.func, second.step$par)
    vcov = solve(t(hss)%*%weight%*%(hss))
    return(list("fs"=first.step, "ss"=second.step, 'vcov'=vcov, 'N'= N))
  } else return(list("fs"=first.step, "ss"=second.step))
}