source("../methods/estimation.R")
source("../methods/selection.R")
library(parallel)

log.lkl <- function(par, y) sum(log(density.function(y,par)))

for(N in SampleSize)
{
  ss=N
  results_list = list()
  Ltest = (length(true.par)):(N)
  
  Ltest = Ltest[Ltest<=max.L]
  
  set.seed(123)
  
  
  modelo = readRDS(paste(mc.name, "_N",ss, ".RDS",sep=""))
  
  #Loads data
  yMat = readRDS( paste("sample_",mc.name,"_N",ss, ".RDS",sep=""))
  
  
  #In this step, we rerun estimators which failed (generated errors) due to convergence issues
  #Checks for caglad 
  val = which(do.call(rbind,lapply(modelo, function(x){sapply(x$caglad_ss, function(y){length(y)!=6})})),arr.ind=T)
  if(nrow(val)>0)
    for(j in 1:nrow(val))
      modelo[[val[j,1]]]$caglad_ss[[val[j,2]]] = lmoment.est(yMat[,val[j,1]], modelo[[val[j,1]]]$caglad_fs[[1]]$fs$par,val[j,2] + length(true.par)-1 , lmoment.analytic = lmoment.analytic, quantile.func = quantile.function,
                                                             density = density.function, grid.length = 2000, par.first.step = modelo[[val[j,1]]]$caglad_fs[[1]]$fs$par, 
                                                             vcov = T, lmoment.deriv.analytic = lmoment.deriv.analytic,
                                                             control = list( "maxit"=500))
  
  #Checks for unbiased
  val = which(do.call(rbind,lapply(modelo, function(x){sapply(x$unbiased_ss, function(y){length(y)!=6})})),arr.ind=T)
  if(nrow(val)>0)
    for(j in 1:nrow(val))
      modelo[[val[j,1]]]$unbiased_ss[[val[j,2]]] = lmoment.est(yMat[,val[j,1]], modelo[[val[j,1]]]$unbiased_fs[[1]]$fs$par,val[j,2] + length(true.par)-1 , lmoment.analytic = lmoment.analytic, quantile.func = quantile.function,
                                                               lmoment.est = "unbiased", 
                                                               vcov=T, lmoment.deriv.analytic = lmoment.deriv.analytic,
                                                               density = density.function, grid.length = 2000, par.first.step = modelo[[val[j,1]]]$unbiased_fs[[1]]$fs$par, control = list( "maxit"=500))
  
  
  #Rerunning some estimators that diverged
  #In this step, we reestimate estimators that diverged (outside the range where L-moments are defined),
  # in the numerical optimization
  
  # Caglad FS estimators
  values = (lapply(modelo, function(x)
    which(!sapply(x$caglad_fs, function(x) x$fs$par[length(x$fs$par)]>=-1))
  ))
  
  for(i in 1:length(values))
    if(length(values[[i]])>0)
      for(vv in values[[i]])
        modelo[[i]]$caglad_fs[[vv]] = lmoment.est(yMat[,i], modelo[[i]]$caglad_fs[[1]]$fs$par, vv+length(true.par)-1, lmoment.analytic = lmoment.analytic, quantile.func = quantile.function, density = density.function, 
                                                  vcov=T, lmoment.deriv.analytic = lmoment.deriv.analytic, par.first.step = modelo[[i]]$caglad_fs[[1]]$fs$par,
                                                  lmoment.est = "caglad",
                                                  weight.matrix = "first",
                                                  grid.length = 2000, control = list( "maxit"=500))
  
  
  
  
  #Unbiased FS
  values = (lapply(modelo, function(x)
    which(!sapply(x$unbiased_fs, function(x) x$fs$par[length(x$fs$par)]>=-1))
  ))
  
  for(i in 1:length(values))
    if(length(values[[i]])>0)
      for(vv in values[[i]])
        modelo[[i]]$unbiased_fs[[vv]] = lmoment.est(yMat[,i], modelo[[i]]$unbiased_fs[[1]]$fs$par, vv+length(true.par)-1, lmoment.analytic = lmoment.analytic, quantile.func = quantile.function, density = density.function, 
                                                    vcov=T, lmoment.deriv.analytic = lmoment.deriv.analytic, par.first.step = modelo[[i]]$unbiased_fs[[1]]$fs$par,
                                                    lmoment.est = "unbiased",
                                                    weight.matrix = "first",
                                                    grid.length = 2000, control = list( "maxit"=500))
  
  
  #Caglad TS
  values = (lapply(modelo, function(x)
    which(!sapply(x$caglad_ss, function(x) x$ss$par[length(x$ss$par)]>=-1))
  ))
  for(i in 1:length(values))
    if(length(values[[i]])>0)
      for(vv in values[[i]])
        modelo[[i]]$caglad_ss[[vv]] = lmoment.est(yMat[,i], modelo[[i]]$caglad_fs[[1]]$fs$par, vv+length(true.par)-1, lmoment.analytic = lmoment.analytic, quantile.func = quantile.function, density = density.function, 
                                                  vcov=T, lmoment.deriv.analytic = lmoment.deriv.analytic, par.first.step = modelo[[i]]$caglad_fs[[1]]$fs$par,
                                                  lmoment.est = "caglad", 
                                                  grid.length = 2000, control = list( "maxit"=500))
  
  
  
  
  
  #Unbiased TS
  values = (lapply(modelo, function(x)
    which(!sapply(x$unbiased_ss, function(x) x$ss$par[length(x$ss$par)]>=-1))
  ))
  for(i in 1:length(values))
    if(length(values[[i]])>0)
      for(vv in values[[i]])
        modelo[[i]]$unbiased_ss[[vv]] = lmoment.est(yMat[,i], modelo[[i]]$unbiased_fs[[1]]$fs$par, vv+length(true.par)-1, lmoment.analytic = lmoment.analytic, quantile.func = quantile.function, density = density.function, 
                                                    vcov=T, lmoment.deriv.analytic = lmoment.deriv.analytic, par.first.step = modelo[[i]]$unbiased_fs[[1]]$fs$par,
                                                    lmoment.est = "unbiased", 
                                                    grid.length = 2000, control = list( "maxit"=500))
  
  
  #Check if something still diverges
  values = (sapply(modelo, function(x) which(!c((x$mle$par[length(x$mle$par)]>=-1),all(sapply(x$caglad_fs, function(x) x$fs$par[length(x$fs$par)]>=-1)),all(sapply(x$caglad_ss, function(x) x$ss$par[length(x$ss$par)]>=-1)),all(sapply(x$unbiased_fs, function(x) x$fs$par[length(x$fs$par)]>=-1)),all(sapply(x$unbiased_ss, function(x) x$ss$par[length(x$ss$par)]>=-1))))))
  
  if(length(unlist(values))>0)
    print("There are still some estimators that diverged") else print("Everything ok!")
  
  
  #Bootstrap
  Bdraws =1000
  Z_draw = matrix(rnorm(Bdraws*(length(true.par))), ncol = Bdraws)
  
  
  grid.length.weight = 2000
  midder = (0:(grid.length.weight-1) + 1:(grid.length.weight))/(2*grid.length.weight)
  mat_middle = outer(midder, rep(1,grid.length.weight))
  
  mat_middle = pmin(mat_middle, t(mat_middle)) - mat_middle*t(mat_middle)
  
  
  power_mat = sapply(0:(max(Ltest)-1), function(l){midder^l})
  
  for(j in 1:Nreps)
  {
    
    
    
    print(j)
    yData = yMat[,j]
    
    xcase = modelo[[j]]$caglad_ss
    
    
    #First-step 

    
    vwd = 1/density.function(quantile.function(midder,xcase[[1]]$fs$par), xcase[[1]]$fs$par)
    vwd[is.nan(vwd)]=0
    mat_mid = mat_middle*(vwd%*%t(vwd))
    
    weight0 = t(power_mat)%*%mat_mid%*%power_mat/grid.length.weight^2
    qdf.grd= grad.qdf(midder, xcase[[1]]$fs$par)
    
    
    cluster = makeCluster(detectCores())
    adder = parLapply( cluster, 1:length(xcase), function(il,  true.par = true.par,  xcase=xcase,mat_middle=mat_middle,
                                                          midder=midder, lmoment.deriv.analytic = lmoment.deriv.analytic,
                                                          quantile.function = quantile.function, density.function=density.function,
                                                          Z_draw=Z_draw, grad.quantile.function=grad.quantile.function, 
                                                          hessian.quantile.function = hessian.quantile.function, power_mat=power_mat,
                                                          grid.length.weight=grid.length.weight,N=N,tau.seq=tau.seq,weight0=weight0#,list_terms=list_terms,lmoment.hessian.analytic=lmoment.hessian.analytic
                                                          ) {
      library(pracma)
      L = length(true.par)+il-1
      y = xcase[[il]]
      
      weightt = weight0[1:L,1:L]
      weight=pinv(weightt)
      
      vwd = 1/density.function(quantile.function(midder,y$ss$par), y$ss$par)
      vwd[is.nan(vwd)]=0
      mat_mid = mat_middle*(vwd%*%t(vwd))
      
      vcov = t(power_mat[,1:L])%*%mat_mid%*%power_mat[,1:L]/grid.length.weight^2
      
      vss = lmoment.deriv.analytic(y$ss$par,L)
      vss[is.nan(vss)] = 0
      Aminus = tryCatch({solve(t(vss)%*%weight%*%vss)}, error=function(e) {pinv(t(vss)%*%weight%*%vss)})
      
      variance = Aminus%*%t(vss)%*%weight%*%vcov%*%weight%*%vss%*%Aminus/N
      
      
 
      
      
      
      distribution = t(tryCatch({chol(variance)},error=function(e) {chol(variance,pivot=T)}))%*%Z_draw
      
      
      table_add = sapply(tau.seq, function(u){
        
        Jacob =       grad.quantile.function(u, y$ss$par)
        Hess =       hessian.quantile.function(u, y$ss$par)
        
        fs.term = Jacob%*%distribution
        
        Ho.term = Hess%*%distribution
        
        
        variance_term = as.numeric((as.numeric(Jacob%*%Aminus%*%t(Jacob)) + 2*t(Ho.term)%*%Aminus%*%t(Jacob))/N) 
        
        t_term = fs.term/sqrt(variance_term)
        
        
        n_ss = fs.term + t(diag(t(distribution)%*%Hess%*%distribution))
        
        Jrep = matrix(rep(Jacob,ncol(distribution)), nrow = ncol(Jacob))
        v_ss = diag(t(Jrep+ Ho.term)%*%Aminus%*%(Jrep + Ho.term))/N
        t_ss = n_ss/sqrt(v_ss)
        
        return(c("mean_tss"=mean(t_ss,na.rm=T), "sd_tss"=sd(t_ss,na.rm=T), 
                 "quantiles_tss"=quantile(t_ss, c(0.025,0.975),na.rm=T),
                 "mean_vss" = mean(v_ss),
                 "sd_vss" = sd(v_ss),
                 "mean_tterm"= mean(t_term,na.rm = T),
                 "sd_tterm" = sd(t_term,na.rm = T),
                 "quantiles_tterm" = quantile(t_term, c(0.025,0.975),na.rm = T),
                 "mean_vterm" = mean(variance_term),
                 "sd_vterm" = sd(variance_term),
                 "sd_orig" = sqrt(Jacob%*%Aminus%*%t(Jacob)/N),
                 "sd_refine" =sqrt(Jacob%*%variance%*%t(Jacob))
                 ))
        
        
      })
      
    }, true.par = true.par,  xcase=xcase,mat_middle=mat_middle,
    midder=midder, lmoment.deriv.analytic = lmoment.deriv.analytic,
    quantile.function = quantile.function, density.function=density.function,
    Z_draw=Z_draw, grad.quantile.function=grad.quantile.function, 
    hessian.quantile.function = hessian.quantile.function, power_mat=power_mat, 
    grid.length.weight=grid.length.weight,N=N, tau.seq=tau.seq,weight0=weight0)
    
    stopCluster(cluster)
    
    results_list[[j]] = adder
    
    if(!dir.exists("bsdraws"))
      dir.create("bsdraws")
    
    
    if(j%%50==0)
      saveRDS(results_list, file = paste("bsdraws/higher",mc.name,"_N",N,"pt",j%/%500+1, ".RDS",sep=""))
    
    
    if(j%%500==0){
      saveRDS(results_list, file = paste("bsdraws/higher",mc.name,"_N",N,"pt",j%/%500, ".RDS",sep=""))
      
      results_list = list()
    }
    
  }
  
  #Deleting a redundant file that the code adds in lines 243-246 when j==Nreps
  file.remove(paste("bsdraws/higher",mc.name,"_N",N,"pt",Nreps%/%500, ".RDS",sep=""))
}
