library(np)
library(pracma)
library(parallel)

library(autodiffr)
ad_setup(JULIA_HOME=julia_path)
library(JuliaCall)
julia_library("SpecialFunctions")

library(quadprog)

lmoment.select <- function(y, par, Lmax, orthogonal = F, uvalues = NULL, lmoment.analytic = NULL, quantile.func=NULL, density.function = NULL,  
                        lmoment.est = "caglad", grid.length = length(y), Nsim = 500, mc.cores = detectCores(),...){
  
  if(orthogonal)
    coefs =  sapply(0:(Lmax-1),function(x){ c((-1)^(x-0:x)*sapply(0:x, function(y){exp(log(choose(x,y))+log(choose(x+y,y)) - log(x+1)*(x+1))}),
                                           rep(0, Lmax-1-x))})
  
  y = y[order(y)]
  N = length(y)
  
  if(lmoment.est=="caglad")
    values = t(sapply(0:(Lmax-1), function(x){
      vals = (((0:N)/N)^(x+1))/(x+1)
      sum(diff(vals)*y)
    })) else if(lmoment.est == "unbiased") values = y%*%sapply(0:(Lmax-1), function(k){sapply(1:N, function(i){exp(log(choose(i-1,k))-log(choose(N,k+1)))})/(k+1)})
  
  
  if(orthogonal)
    l.est = values%*%coefs else l.est = values
  
  if(is.null(lmoment.analytic))
  {
    
    midpointer = (0:(grid.length-1) + 1:(grid.length))/(2*grid.length)
    
    powers = sapply(0:(Lmax-1), function(x){
      midpointer^x})
    
    lmoment.func <- function(theta){
      vvv = colMeans(sapply(midpointer, function(u){quantile.func(u,theta)})*powers)
      if(orthogonal)
        return(vvv%*%coefs) else return(vvv)
    }
  } else {
    lmoment.func <- function(theta){
      lmoment.analytic(theta, Lmax)
    }
  }
  Lmin = length(par)
  objective <- function(theta){t(lmoment.func(theta)[1:Lmin]-l.est[1:Lmin])%*%(lmoment.func(theta)[1:Lmin]  -l.est[1:Lmin])}

  first.step = optim(par, objective, method = "BFGS", ...) 
  
  l.theta.fs = lmoment.func(first.step$par)
  
  l.est.mat = do.call(cbind,lapply(1:Nsim, function(j)
  {
    ysim = quantile.func(runif(N), first.step$par)
      
    
    ysim = ysim[order(ysim)]
    
    if(lmoment.est=="caglad")
      values_sim = t(sapply(0:(Lmax-1), function(x){
        vals = (((0:N)/N)^(x+1))/(x+1)
        sum(diff(vals)*ysim)
      })) else if(lmoment.est == "unbiased") values_sim = ysim%*%sapply(0:(Lmax-1), function(k){sapply(1:N, function(i){exp(log(choose(i-1,k))-log(choose(N,k+1)))})/(k+1)})
      
    if(orthogonal)
      l.est.sim = values_sim%*%coefs else l.est.sim = values_sim   
    
    t(l.est.sim)
    
  
  } ))
  
  
  h1 = - ad_jacobian(lmoment.func, first.step$par)
  
  h2 = lapply(1:Lmax, function(j)
    {
    ff <- function(par) lmoment.func(par)[j]
    
    - ad_hessian(ff, first.step$par)
  })
  
  #We start with the first step  simulated representation
  fs.rep =  - solve(t(h1[1:Lmin,])%*%h1[1:Lmin,])%*%t(h1[1:Lmin,])%*%(l.est.mat[1:Lmin,]- l.theta.fs[1:Lmin])
  
  #Now calculate optimal weights using the parametric estimator
  grid.length.weight = grid.length
  midder = (0:(grid.length.weight-1) + 1:(grid.length.weight))/(2*grid.length.weight)
  
  vwd = sapply(midder, function(u){1/density.function(quantile.func(u,first.step$par), first.step$par)})
  
  mat_middle = outer(midder, rep(1,grid.length.weight))
  
  mat_middle = pmin(mat_middle, t(mat_middle)) - mat_middle*t(mat_middle)
  
  power_mat = sapply(0:(Lmax-1), function(l){midder^l})
  
  weight = t(power_mat)%*%(mat_middle*((vwd%*%t(vwd))))%*%power_mat/grid.length.weight^2
  
  if(orthogonal)
    weight = t(coefs)%*%weight%*%coefs
  
  
  #Calculating derivatives of optimal weighting matrix
  jac_weight = ad_jacobian(function(par) 1/density.function(quantile.func(midder,par), par), first.step$par) 
  
  mat_fs_deriv = lapply(1:ncol(jac_weight), function(j){
    jjj = t(power_mat)%*%(mat_middle*(vwd%*%t(jac_weight[,j]) + jac_weight[,j]%*%t(vwd)))%*%power_mat/grid.length.weight^2 
    
    if(orthogonal)
      jjj = t(coefs)%*%jjj%*%coefs
    
    return(jjj)
  })
  
  Lgrid = seq(Lmin, Lmax, 1)
  
  cluster = makeCluster(mc.cores)
  Valores_Theta = parLapply(cl = cluster,  X = Lgrid, fun = function(L, weight, l.est.mat, l.theta.fs, h1, h2, par, mat_fs_deriv, fs.rep)
  {
    library(pracma)
    #print(L)
    Sigma = pinv(weight[1:L,1:L])
    Sigma.star = pinv(t(h1[1:L,])%*%Sigma%*%h1[1:L,])
    
    diff.H = -(l.est.mat[1:L,]- l.theta.fs[1:L])
    
    Mo.inv = cbind(Sigma.star, - Sigma.star%*%t(h1[1:L,])%*%Sigma)
    Mo.inv = rbind(Mo.inv,
                   cbind(-Sigma%*%h1[1:L,]%*%Sigma.star, - Sigma + Sigma%*%h1[1:L,]%*%Sigma.star%*%t(h1[1:L,])%*%Sigma))
    
    m0 = rbind( matrix(0, nrow = length(par), ncol = ncol(diff.H)),diff.H[1:L,])
    Theta1 = - Mo.inv%*%m0
    
    
    DeltaM = lapply(1:(length(par)+L), function(j)
    {
      if(j<=length(par))
      {
        hff = (sapply(h2[1:L], function(hf) hf[,j]))
        
        matter = cbind(matrix(0,length(par),length(par)), -hff)
        matter = rbind(matter, cbind(-t(hff), matrix(0,L,L)))
      } else matter = rbind(cbind(-h2[[j-length(par)]],matrix(0, nrow = length(par), ncol  = L)), matrix(0, nrow = L, ncol = L+length(par)))
      matter
    })
    
    Theta2 = sapply(1:Nsim, function(j){
      
      if(L> length(par))
      {
        M.diff = Reduce("+", lapply(1:length(par), function(k) - mat_fs_deriv[[k]][1:L, 1:L]*fs.rep[k,j]))
      
        pp = c(rep(0, length(par)), M.diff%*%Theta1[(length(par)+1):nrow(Theta1),j])
        
        part2 =  -Mo.inv%*%pp  - (Mo.inv/2)%*%Reduce("+",lapply(1:(L+length(par)), function(s) Theta1[s, j]*DeltaM[[s]]%*%Theta1[,j]))
        
      } else  part2 =  - (Mo.inv/2)%*%Reduce("+",lapply(1:(L+length(par)), function(s) Theta1[s, j]*DeltaM[[s]]%*%Theta1[,j]))
    })
    
    Thetaf = Theta1+Theta2
    
    return(Thetaf[1:length(par),])
  }, weight = weight, l.est.mat = l.est.mat, l.theta.fs = l.theta.fs, h1=h1, h2 = h2, par = par, mat_fs_deriv = mat_fs_deriv, fs.rep = fs.rep)
  
  
  stopCluster(cluster)
  
  if(is.null(uvalues))
    return(list("VT" = Valores_Theta, "first.step" = first.step)) else {
      
      Lvals = Lgrid[sapply(uvalues, function(u){
        qff = function(par) quantile.func(u, par)
        
        qd = ad_grad(qff, first.step$par)
        qh = ad_hessian(qff, first.step$par)
        
        valores = sapply(Valores_Theta, function(Theta) qd%*%Theta + diag(t(Theta)%*%qh%*%Theta))
        
        which.min(colMeans(valores^2))
      })]
      
      return(list("Lvals" = Lvals, "first.step" = first.step))
    }
  
  
 
}


lmoment.lasso <- function(y, par, Lmax, orthogonal = F, lmoment.analytic = NULL, quantile.func=NULL, density.function = NULL,  
                           lmoment.est = "caglad",   weight.matrix = "par", grid.length = length(y), max.iter = 50, tol.iter = 0.01, step.iter = 0.01, ...){
  
  if(orthogonal)
    coefs =  sapply(0:(Lmax-1),function(x){ c((-1)^(x-0:x)*sapply(0:x, function(y){exp(log(choose(x,y))+log(choose(x+y,y)) - log(x+1)*(x+1))}),
                                              rep(0, Lmax-1-x))})
  
  y = y[order(y)]
  N = length(y)
  
  if(lmoment.est=="caglad")
    values = t(sapply(0:(Lmax-1), function(x){
      vals = (((0:N)/N)^(x+1))/(x+1)
      sum(diff(vals)*y)
    })) else if(lmoment.est == "unbiased") values = y%*%sapply(0:(Lmax-1), function(k){sapply(1:N, function(i){exp(log(choose(i-1,k))-log(choose(N,k+1)))})/(k+1)})
  
  
  if(orthogonal)
    l.est = values%*%coefs else l.est = values
  
  if(is.null(lmoment.analytic))
  {
    
    midpointer = (0:(grid.length-1) + 1:(grid.length))/(2*grid.length)
    
    powers = sapply(0:(Lmax-1), function(x){
      midpointer^x})
    
    lmoment.func <- function(theta){
      vvv = colMeans(sapply(midpointer, function(u){quantile.func(u,theta)})*powers)
      if(orthogonal)
        return(vvv%*%coefs) else return(vvv)
    }
  } else {
    lmoment.func <- function(theta){
      lmoment.analytic(theta, Lmax)
    }
  }
  Lmin = length(par)
  objective <- function(theta){t(lmoment.func(theta)[1:Lmin]-l.est[1:Lmin])%*%(lmoment.func(theta)[1:Lmin]  -l.est[1:Lmin])}
  
  first.step = optim(par, objective, method = "BFGS", ...) 
  
  l.theta.fs = lmoment.func(first.step$par)
  
  
  h1 = - ad_jacobian(lmoment.func, first.step$par)
  
  h2 = lapply(1:Lmax, function(j)
  {
    ff <- function(par) lmoment.func(par)[j]
    
    - ad_hessian(ff, first.step$par)
  })

  #Now calculate optimal weights using the parametric estimator
  
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
  
  power_mat = sapply(0:(Lmax-1), function(l){midder^l})
  
  weight = t(power_mat)%*%(mat_middle*((vwd%*%t(vwd))))%*%power_mat/grid.length.weight^2
  
  if(orthogonal)
    weight = t(coefs)%*%weight%*%coefs
  
  t0 = 1.1*sqrt(N*qnorm(1- 0.1/(log(max(Lmax), N)*4*Lmax*length(par))))/N
  

  A_prog = rbind(cbind(weight, - weight),cbind(-weight, weight))
  
  A_prog = nearest_spd(A_prog)

  
  #Calculating derivatives of optimal weighting matrix
  jac_weight = ad_jacobian(function(par) 1/density.function(quantile.func(midder,par), par), first.step$par) 
  
  mat_fs_deriv = lapply(1:ncol(jac_weight), function(j){
    jjj = t(power_mat)%*%(mat_middle*(vwd%*%t(jac_weight[,j]) + jac_weight[,j]%*%t(vwd)))%*%power_mat/grid.length.weight^2 
    
    if(orthogonal)
      jjj = t(coefs)%*%jjj%*%coefs
    
    return(jjj)
  })
  
  bread.fs =  solve(t(h1[1:Lmin,])%*%h1[1:Lmin,])%*%t(h1[1:Lmin,])
  meat.fs = t(power_mat[,1:Lmin])%*%(mat_middle*((vwd%*%t(vwd))))%*%power_mat[,1:Lmin]/grid.length.weight^2

  fs.rep =  bread.fs%*%meat.fs%*%t(bread.fs)
  
  
  upA = sapply(1:Lmax, function(j){
    
    vvK = sapply(mat_fs_deriv, function(jk) jk[j,])
    
    sqrt(max(diag(vvK%*%fs.rep%*%t(vvK)))) 
  })
  
  upA[1:length(par)] = 0
  
  Matt = sapply(1:length(par), function(p){
    
    bb = h1[,p]
    

    upB = sapply(1:Lmax, function(j){
      vvK = h2[[j]][p,]
      
      sqrt(t(vvK)%*%fs.rep%*%vvK )
    })
    
    
    upB[1:length(par)] = 0
    
    penalty_low =  sum(t0*(upA*0 + upB))
    
    penalty_upper = sum(t0*(upA*Lmax + upB))
    
    f_opt <- function(wL) (wL - sum(t0*(upA*wL + upB)))^2
    
    opter = optimize(f_opt,lower = penalty_low, upper = penalty_upper)
    
    
    
    bb_prog = c(bb, - bb) 
    
    penalty_base = t0*(upA*opter$minimum + upB)
    
    
    qpp = tryCatch({solve.QP(Dmat = A_prog, dvec = +bb_prog - rep(penalty_base, 2), Amat = diag(nrow(A_prog)))},
                    error = function(e)
                    {
                      Dmat = A_prog + 1e-10*diag(nrow(A_prog))
                      #print(Dmat)
                      solve.QP(Dmat = Dmat , dvec = +bb_prog - rep(penalty_base, 2), Amat = diag(nrow(A_prog)))
                    })

    lambda_new = qpp$solution

    lambda_new[qpp$iact] = 0
    lambda_new = lambda_new[1:nrow(weight)] - lambda_new[(nrow(weight)+1):(2*nrow(weight))]

    lambda = lambda_new

    # itt = 0
    # 
    # repeat{
    # 
    #   qpp = solve.QP(Dmat = A_prog, dvec = +bb_prog - rep(penalty_base, 2), Amat = diag(nrow(A_prog)))
    # 
    #   lambda_new = qpp$solution
    # 
    #   lambda_new[qpp$iact] = 0
    #   lambda_new = lambda_new[1:nrow(weight)] - lambda_new[(nrow(weight)+1):(2*nrow(weight))]
    # 
    # 
    # 
    #   if(itt==0)
    #   {
    #     lambda_old = 0
    #     lambda = lambda_new
    #   } else {
    #     lambda_old = lambda
    #     lambda = (1-step.iter)*lambda + step.iter*lambda_new
    #   }
    # 
    #   if(itt==0)
    #     err = tol.iter*2 else err =  mean(abs(lambda-lambda_old))
    # 
    #   print(itt)
    #   print(err)
    #   itt = itt+1
    # 
    #   if(err<tol.iter|itt> max.iter)
    #     break
    # 
    # 
    #   upAlt = sapply(1:Lmax, function(j){
    # 
    #     vvK = sapply(mat_fs_deriv, function(jk) jk[j,])
    # 
    #     sqrt(t(lambda)%*%vvK%*%fs.rep%*%t(vvK)%*%lambda)
    #   })
    # 
    #   penalty_base = t0*(upAlt + upB)
    # }

    return(lambda)
  })
  
  objective <- function(theta){((lmoment.func(theta)-l.est)%*%(Matt))%*%t((lmoment.func(theta)-l.est)%*%(Matt))}
  
  second.step = optim(par, objective, method = "BFGS", ...) 
  
  
  keeper = rowSums(abs(Matt))>0
  
  if(sum(keeper)>0)
  {
    if(sum(keeper)==1)
      weight_post = as.matrix(1/weight[keeper,keeper]) else weight_post = pinv(weight[keeper,keeper])
      
    objective <- function(theta){t((lmoment.func(theta)-l.est)[keeper])%*%weight_post%*%((lmoment.func(theta)-l.est)[keeper])}
    
    post.second.step =    optim(par, objective, method = "BFGS",...) 
  } else post.second.step = first.step
  
  return(list("fs" = first.step, "ss" = second.step, "ps" = post.second.step, "kept" = sum(keeper)))
}
  