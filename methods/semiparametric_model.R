#This code estimates a parametric model for the distribution of the error term epsilon of the semiparametric model:
# Y = a(X) + b(X)W + epsilon,
#where a and b are identified by E[epsilon|X]=0 and E[S epsilon|X]=0, with S an instrumental variable.
# The functions a and b are estimated with the linear instrumental forests algorithm of Athey et al
# The parametric model for epsilon is estimated through the estimator in Supplemental Appendix J of Alvarez et al
library(grf)

lmoment.semiparametric <- function(train, holdout, y_pos, x_pos, w_pos, s_pos, par, L, orthogonal = F, grid_inv=nrow(holdout), lmoment.analytic = NULL, quantile.func=NULL, density.function = NULL,  
                           lmoment.est = "caglad", grid.length = nrow(holdout), two.step = T,
                           vcov.est = "par",bounds=NULL,model = NULL,...){
  
  if(is.null(model))
  {
  models <- lapply(c("train", "holdout"),
                   function(j) {
                     data = get(j)
                     slope = instrumental_forest(as.matrix(data[,x_pos]), as.matrix(data[,y_pos]), as.matrix(data[,w_pos]), as.matrix(data[,s_pos]))
                     intercept = regression_forest(as.matrix(data[,x_pos]), as.numeric(data[,y_pos]) - predict(slope)[,1]*as.numeric(data[,w_pos]))
                      
                    return(list("slope"=slope,"intercept"=intercept))
                     
                   }
                     )
  names(models) = c("train", "holdout")
  
  predict.mod <- function(model, j){
    dataset = get(j)
    
    if(model==j)
      return(predict(models[[model]]$intercept)[,1] + (predict(models[[model]]$slope))[,1]*dataset[,w_pos]) else return(predict(models[[model]]$intercept, newdata = dataset[,x_pos])[,1] + (predict(models[[model]]$slope, newdata = dataset[,x_pos])[,1])*dataset[,w_pos])
  }
  
  res_oos = holdout[,y_pos] - predict.mod("train", "holdout")
  
  res_iis = holdout[,y_pos] - predict.mod("holdout", "holdout")
  
  res_rev =  train[,y_pos] - predict.mod("holdout", "train")
  
  bX = lm_forest(as.matrix(holdout[,x_pos]), as.matrix(holdout[,s_pos]), as.matrix(holdout[,w_pos]))
  aX =  regression_forest(as.matrix(holdout[,x_pos]), as.numeric(holdout[,s_pos] - (predict(bX)$predictions[,1,1])*holdout[,w_pos])) 
  
  S_rank = (train[,s_pos] - predict(aX, newdata = train[x_pos])[,1])/predict(bX, newdata = as.matrix(train[,x_pos]))$predictions[,1,1]
  
  grid_u = seq(1/grid_inv, 1 - 1/grid_inv, by = 1/grid_inv)
  
  qe = quantile(res_iis, grid_u)
  
  upper_q = quantile(res_iis, grid_u + 1/(2*grid_inv))
  lower_q = quantile(res_iis, grid_u - 1/(2*grid_inv))
  adjust = (upper_q - lower_q)/(1/(grid_inv))
  
  list_results = lapply(1:length(grid_u), function(rr)
  {
    print(rr)
    up = upper_q[rr]
    down = lower_q[rr]
    
    slope_upper = lm_forest(X = as.matrix(holdout[,x_pos]), Y = 1*(res_iis <= up),
                            W = as.matrix(holdout[,w_pos]))
    slope_lower = lm_forest(X = as.matrix(holdout[,x_pos]), Y = 1*(res_iis <= down),
                            W = as.matrix(holdout[,w_pos]))
    
    intercept_upper = regression_forest(as.matrix(holdout[,x_pos]), 1*(res_iis <= up) - predict(slope_upper)$predictions[,1,1]*holdout[,w_pos])
    intercept_lower = regression_forest(as.matrix(holdout[,x_pos]), 1*(res_iis <= down) - predict(slope_upper)$predictions[,1,1]*holdout[,w_pos])
    
    rs = (predict(intercept_upper, newdata = as.matrix(train[,x_pos]))[,1] - predict(intercept_lower, newdata = as.matrix(train[,x_pos]))[,1])/(up-down) +     S_rank*(predict(slope_upper, newdata = as.matrix(train[,x_pos]))$predictions[,1,1] - predict(slope_lower, newdata = as.matrix(train[,x_pos]))$predictions[,1,1])/(up-down)
    
    return(rs)
}
    
    )
  
  
  mat = do.call(rbind, list_results)
  
  mat = adjust*mat*(rep(1,length(grid_u))%*%t(res_rev))
  
  corr = rowMeans(mat)
  
  } else {
    res_oos = model$res_oos
    corr = model$corr
  }
  

  
  if(orthogonal)
    coefs =  sapply(0:(L-1),function(x){ c((-1)^(x-0:x)*sapply(0:x, function(y){exp(log(choose(x,y))+log(choose(x+y,y)) - log(x+1)*(x+1))}),
                                           rep(0, L-1-x))})
  
  
  val = lapply(list(res_oos, corr), function(y) {

  y = y[order(y)]
  N = length(y)
  
  if(lmoment.est=="caglad")
    values = t(sapply(0:(L-1), function(x){
      vals = (((0:N)/N)^(x+1))/(x+1)
      sum(diff(vals)*y)
    })) else if(lmoment.est == "unbiased") values = y%*%sapply(0:(L-1), function(k){sapply(1:N, function(i){exp(log(choose(i-1,k))-log(choose(N,k+1)))})/(k+1)})
  
  
  if(orthogonal)
    l.est = values%*%coefs else l.est = values
  
  return(l.est)
  })
  
  l.est = val[[1]]
  
  A_add  = val[[2]]
  
  y = res_oos
  
  y = y[order(y)]
  
  N = length(y)
  

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
  
  objective <- function(theta){
    (lmoment.func(theta)-l.est -A_add)%*%t(lmoment.func(theta) -l.est-A_add)}
  
  if(is.null(bounds))
    first.step = optim(par, objective, method = "BFGS", ...) else first.step = optim(par, objective, method = "L-BFGS-B", lower = bounds$lower, upper = bounds$upper, ...)
  
  #print(first.step)
  if(vcov.est == "semi")
  {
    grid.length.vcov = N
    midder = (0:(grid.length.vcov-1) + 1:(grid.length.vcov))/(2*grid.length.vcov)
    vwd = 1/sapply(y, function(x){density.function(x, first.step$par)})
    
  } else if(vcov.est == "par")
  {
    grid.length.vcov = grid.length
    midder = (0:(grid.length.vcov-1) + 1:(grid.length.vcov))/(2*grid.length.vcov)
    vwd = sapply(midder, function(u){1/density.function(quantile.func(u,first.step$par), first.step$par)})
    
  } else if(vcov.est == "np")
  {
    grid.length.vcov = N
    midder = (0:(grid.length.vcov-1) + 1:(grid.length.vcov))/(2*grid.length.vcov)
    
    bw = npudensbw(y)
    vwd = 1/fitted(npudens(bw))
  }
  
  
  mat_middle = outer(midder, rep(1,grid.length.vcov))
  
  mat_middle = pmin(mat_middle, t(mat_middle)) - mat_middle*t(mat_middle)
  
  power_mat = sapply(0:(L-1), function(l){midder^l})
  
  weight = t(power_mat)%*%(mat_middle*((vwd%*%t(vwd))))%*%power_mat/grid.length.vcov^2
  
  if(orthogonal)
    weight = t(coefs)%*%weight%*%coefs
  
  middle_sand = weight
  
  #print(vwd)
  if(!two.step)
  {
    second.step = first.step
    
    weight = diag(L)
    
  } else {
    weight = pinv(middle_sand)
    objective <- function(theta){(lmoment.func(theta)-l.est -A_add)%*%weight%*%(t(lmoment.func(theta)  -l.est -A_add))}
    
    #print(first.step)
    if(is.null(bounds))
      second.step = tryCatch({optim(first.step$par, objective, method = "BFGS", ...)},
                             error = function(e){optim(par, objective, method = "BFGS", ...)}) else second.step = optim(first.step$par, objective, method = "L-BFGS-B", lower = bounds$lower, upper = bounds$upper, ...)
    
  }
  
  h1 = - ad_jacobian(lmoment.func, second.step$par)
  variance_estimator = solve(t(h1)%*%weight%*%h1)%*%(t(h1)%*%weight%*%middle_sand%*%weight%*%h1)%*%solve(t(h1)%*%weight%*%h1)

  if(is.null(model))
    model.add = list("res_oos" = res_oos, "res_iis" = res_iis, "res_rev" = res_rev, 'corr' = corr, 'mat' = mat,
                                       'jacob' = h1, 'weight'=weight, 'sand'= middle_sand,
                                       'l.est'=l.est, 'A_add' = A_add, 'fs'=first.step) else  model.add = list("res_oos" = res_oos,
                        'jacob' = h1, 'weight'=weight, 'sand'= middle_sand,
                        'l.est'=l.est, 'A_add' = A_add, 'fs'=first.step)
                   
    
  return(list("est" = second.step, "vcov" = variance_estimator,"N"=N,'model'=model.add))
}