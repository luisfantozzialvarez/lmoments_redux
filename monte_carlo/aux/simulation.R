source("../methods/estimation.R")
source("../methods/mle_tilting.R")
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
    
    mle = optim(true.par, log.lkl, method = "BFGS", control =  list("fnscale"=-1, "maxit"=500), y = yData, hessian = T)
    
    #Trimmed MLE
    trimmed_mle <- lapply(tau.seq[tau.seq>0.5], function(x)
      {
      y_trim = yData[yData < quantile.function(x,mle$par)]
      
      
      tryCatch({optim(true.par, log.lkl, method = "BFGS", control =  list("fnscale"=-1, "maxit"=500), y = y_trim)},error=function(e) NULL)
      
    })
    
    
    #Tilted_MLE
    tilted_mle <- lapply(tau.seq[tau.seq>0.5], function(x)
    {
      tryCatch({empirically.tilted(yData, density.function, true.par, delta = -log(x))}, error=function(e) NULL)
    })
    
    
    caglad_fs =  mclapply(Ltest, function(L){lmoment.est(yData, true.par, L, lmoment.analytic = lmoment.analytic, quantile.func = quantile.function,
                                                         density = density.function, grid.length = 2000, weight.matrix = "first", control = list( "maxit"=500))},
                          mc.cores = detectCores())
    
    caglad_ss = mclapply(Ltest, function(L){lmoment.est(yData, true.par, L, lmoment.analytic = lmoment.analytic, quantile.func = quantile.function,
                density = density.function, grid.length = 2000, par.first.step = caglad_fs[[1]]$fs$par, vcov=T, lmoment.deriv.analytic = lmoment.deriv.analytic, control = list( "maxit"=500))},
                mc.cores = detectCores())
    
    unbiased_fs = mclapply(Ltest, function(L){lmoment.est(yData, true.par, L, lmoment.analytic = lmoment.analytic, quantile.func = quantile.function,
                                                     density = density.function,
                                                     lmoment.est = "unbiased", grid.length = 2000,  weight.matrix = "first", control = list( "maxit"=500))},
                      mc.cores = detectCores())
    
    unbiased_ss = mclapply(Ltest, function(L){lmoment.est(yData, true.par, L, lmoment.analytic = lmoment.analytic, quantile.func = quantile.function,
                                                          density = density.function,
                                                          lmoment.est = "unbiased", grid.length = 2000, par.first.step = unbiased_fs[[1]]$fs$par, vcov=T, lmoment.deriv.analytic = lmoment.deriv.analytic, control = list( "maxit"=500))},
                           mc.cores = detectCores())
   
    results_list[[j]] = list("mle"=mle,"caglad_fs" = caglad_fs, "caglad_ss" = caglad_ss, 
                             "unbiased_fs" = unbiased_fs, "unbiased_ss" = unbiased_ss,
                             'trimmed_mle'=trimmed_mle, 'tilted_mle'=tilted_mle)
    

    if(j%%50==0)
      saveRDS(results_list, file = paste(mc.name,"_N",N, ".RDS",sep=""))
  
  }
  saveRDS(results_list, file = paste(mc.name,"_N",N, ".RDS",sep=""))

}
  