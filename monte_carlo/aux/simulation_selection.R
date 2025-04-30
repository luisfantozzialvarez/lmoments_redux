source("../methods/estimation.R")
source("../methods/selection.R")
library(parallel)

log.lkl <- function(par, y) sum(log(density.function(y,par)))

for(N in SampleSize)
{
  print(N)
  results_list = list()
  Ltest = (length(true.par)):(N)
  
  Ltest = Ltest[Ltest<=max.L]
  
  yMat = readRDS( paste("sample_",mc.name,"_N",N, ".RDS",sep=""))
  
  set.seed(123)
  for(j in 1:Nreps)
  {
 
    
    
    print(j)
    yData = yMat[,j]
    
    mle = optim(true.par, log.lkl, method = "BFGS", control =  list("fnscale"=-1, "maxit"=500), y = yData)
    
    mat_select_caglad = tryCatch({lmoment.select(yData, true.par, max(Ltest), orthogonal = F, uvalues = tau.seq, lmoment.analytic = lmoment.analytic, quantile.func=quantile.function, density.function = density.function,  
                               lmoment.est = "caglad", grid.length = 2000, Nsim = 1000,
                               lmoment.deriv.analytic =lmoment.deriv.analytic,
                               lmoment.hessian.analytic = lmoment.hessian.analytic, grad.qdf.analytic = grad.qdf, grad.qf.analytic = grad.quantile.function, hessian.qf.analytic=hessian.quantile.function, 
                               control = list( "maxit"=500),mc.cores= detectCores() )},
                               error = function(e){
                                 print(e)
                                 lmoment.select(yData, mle$par, max(Ltest), orthogonal = F, uvalues = tau.seq, lmoment.analytic = lmoment.analytic, quantile.func=quantile.function, density.function = density.function,  
                                                lmoment.est = "caglad", grid.length = 2000, Nsim = 1000,
                                                lmoment.deriv.analytic =lmoment.deriv.analytic,
                                                lmoment.hessian.analytic = lmoment.hessian.analytic, grad.qdf.analytic = grad.qdf, grad.qf.analytic = grad.quantile.function, hessian.qf.analytic=hessian.quantile.function, 
                                                control = list( "maxit"=500),mc.cores= detectCores() )
                                 
                               })
    

    print(mat_select_caglad$Lvals)
    
 
    caglad_ss = mclapply(mat_select_caglad$Lvals, function(L){tryCatch({lmoment.est(yData, true.par, L, lmoment.analytic = lmoment.analytic, quantile.func = quantile.function,
                                                                              density = density.function,
                                                                              lmoment.est = "caglad", grid.length = 2000, par.first.step = mat_select_caglad$first.step$par,
                                                                              control = list( "maxit"=500))},
                                                                        error = function(e) {
                                                                          print(e)
                                                                          lmoment.est(yData, mle$par, L, lmoment.analytic = lmoment.analytic, quantile.func = quantile.function,
                                                                                      density = density.function,
                                                                                      lmoment.est = "caglad", grid.length = 2000, par.first.step = mat_select_caglad$first.step$par,        
                                                                                      control = list( "maxit"=500))
                                                                        })},
                           mc.cores = 1)
    
    
    
    lasso.select = tryCatch({tryCatch({lmoment.lasso(yData, true.par, 2*max(Ltest), orthogonal = F, lmoment.analytic = lmoment.analytic, quantile.func=quantile.function, density.function = density.function,
                              lmoment.est = "caglad",   weight.matrix = "par", grid.length = 2000, max.iter = 10, tol.iter = 0.01, step.iter = 0.01, 
                              lmoment.deriv.analytic =lmoment.deriv.analytic,
                              lmoment.hessian.analytic = lmoment.hessian.analytic, grad.qdf.analytic = grad.qdf, 
                              control = list("maxit" = 500))},
                            error = function(e){
                              lmoment.lasso(yData, mle$par, 2*max(Ltest), orthogonal = F, lmoment.analytic = lmoment.analytic, quantile.func=quantile.function, density.function = density.function,
                                            lmoment.est = "caglad",   weight.matrix = "par", grid.length = 2000, max.iter = 10, tol.iter = 0.01, step.iter = 0.01, 
                                            lmoment.deriv.analytic =lmoment.deriv.analytic,
                                            lmoment.hessian.analytic = lmoment.hessian.analytic, grad.qdf.analytic = grad.qdf, 
                                            control = list("maxit" = 500))
                            })},
                            error=function(e){
                              tryCatch({lmoment.lasso(yData, true.par, max(Ltest), orthogonal = F, lmoment.analytic = lmoment.analytic, quantile.func=quantile.function, density.function = density.function,
                                                      lmoment.est = "caglad",   weight.matrix = "par", grid.length = 2000, max.iter = 10, tol.iter = 0.01, step.iter = 0.01, 
                                                      lmoment.deriv.analytic =lmoment.deriv.analytic,
                                                      lmoment.hessian.analytic = lmoment.hessian.analytic, grad.qdf.analytic = grad.qdf, 
                                                      control = list("maxit" = 500))},
                                       error = function(e){
                                         lmoment.lasso(yData, mle$par, max(Ltest), orthogonal = F, lmoment.analytic = lmoment.analytic, quantile.func=quantile.function, density.function = density.function,
                                                       lmoment.est = "caglad",   weight.matrix = "par", grid.length = 2000, max.iter = 10, tol.iter = 0.01, step.iter = 0.01, 
                                                       lmoment.deriv.analytic =lmoment.deriv.analytic,
                                                       lmoment.hessian.analytic = lmoment.hessian.analytic, grad.qdf.analytic = grad.qdf, 
                                                       control = list("maxit" = 500))
                                       })
                            })

    results_list[[j]] = list("mle"=mle,"caglad_ss" = caglad_ss, "select_caglad" = mat_select_caglad,
                             "lasso_select" = lasso.select)

    if(j%%50==0)
      saveRDS(results_list, file = paste("select",mc.name,"_N",N, ".RDS",sep=""))
    
 
  }
  saveRDS(results_list, file = paste("select",mc.name,"_N",N, ".RDS",sep=""))
}