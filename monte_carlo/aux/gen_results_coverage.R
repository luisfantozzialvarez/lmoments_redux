library(xtable)
library(parallel)
source("../methods/estimation.R")

table = c()

for(ss in SampleSize)
{
  
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
  {print("There are still some estimators that diverged") 
    keeper = !sapply(values, function(x) length(x)>0)
    
    #Dropping diverged draws
    keeper = !sapply(values, function(x) length(x)>0)
    modelo = modelo[keeper]
    
    yMat = yMat[,keeper]
    
    print(paste('Dropped', sum(!keeper), 'simulations that diverged'))
    
  }else print("Everything ok!")
  
  
  if(file.exists(paste("bsdraws/higher",mc.name,"_N",ss,".RDS",sep="")))
    bs_draws = readRDS(paste("bsdraws/higher",mc.name,"_N",ss,".RDS",sep="")) else {
      files= list.files('bsdraws')[grepl(paste("higher",mc.name,"_N",ss,"pt",sep=""),list.files('bsdraws'))]
      files = files[order(files)]
      bs_draws = c()
      for(ff in files)
      {
        aaa = readRDS(paste('bsdraws/',ff,sep=""))
        aaa = aaa[sapply(aaa, function(x) !is.null(x))]
        
        bs_draws= c(bs_draws,aaa)
      }
    }
  
  if(length(unlist(values))>0)
  {
    bs_draws = bs_draws[keeper]
  }
 
  
  
  coverage = 0.95
  

  
  for(tau in tau.seq)
  {
    ff = quantile.function(tau, true.par)
    
    
    bias = sapply(bs_draws, function(y)
      sapply(y, function(x)
      {
        x["mean_tterm",which(tau==tau.seq)]
      }))
    
    dof_adj = sapply(bs_draws, function(y)
      sapply(y, function(x)
      {
        vterm=x["sd_tterm",which(tau==tau.seq)]^2
        ifelse(vterm<=1,Inf, 2*vterm/(vterm-1))
      }))
    
    qlow = sapply(bs_draws, function(y)
      sapply(y, function(x)
      {
        x[grepl("quantiles_tterm",rownames(x)),which(tau==tau.seq)][1]
      }))
    
    qhigh = sapply(bs_draws, function(y)
      sapply(y, function(x)
      {
        x[grepl("quantiles_tterm",rownames(x)),which(tau==tau.seq)][2]
      }))
    

  
    
    
    
    mle_se = sapply(modelo, function(x) {
                    jacob = grad.quantile.function(tau, x$mle$par);
                    sqrt(-(jacob)%*%solve(x$mle$hessian/ss)%*%(t(jacob)))/sqrt(ss)})
    
    ts_se = sapply(modelo, function(x){
      sapply( x$caglad_ss, function(y) { 
        jacob = grad.quantile.function(tau, y$ss$par)
        tryCatch({sqrt((jacob)%*%solve(y$meat.ts)%*%(t(jacob)))/sqrt(ss)},error= function(e)
        {sqrt((jacob)%*%pinv(y$meat.ts)%*%(t(jacob)))/sqrt(ss)})
    }
        )})
    
    

    
   
    mle_est = sapply(modelo, function(x){
      quantile.function(tau,x$mle$par)
    })
    
    ts_est =  sapply(modelo, function(x){sapply(x$caglad_ss, function(y){quantile.function(tau,y$ss$par)})})
    

    
      
    
    
    crit_norm = qnorm(1-(1-coverage)/2)
    
    cov_mle =   mean((ff>=mle_est-crit_norm*mle_se)&(ff<=mle_est+crit_norm*mle_se))
    cov_mle_true =   mean((ff>=mle_est-crit_norm*sd(mle_est))&(ff<=mle_est+crit_norm*sd(mle_est)))
    ci_length_mle_est = median(mle_se)/sd(mle_est)
   
    
    cov_lmoment =  rowMeans((ff>=ts_est-crit_norm*ts_se)&(ff<=ts_est+crit_norm*ts_se),na.rm=T)
    cov_lmoment_true =  rowMeans((ff>=ts_est-crit_norm*apply(ts_est,1,sd,na.rm=T))&(ff<=ts_est+crit_norm*apply(ts_est,1,sd,na.rm=T)),na.rm=T)
    ci_length_true = apply(ts_est,1,sd,na.rm=T)/sd(mle_est)
    ci_length_est = apply(ts_se,1,median,na.rm=T)/sd(mle_est)
  
    
    
    cov_lmoment_qt =  rowMeans((ff>=ts_est -qhigh*ts_se)&(ff<=ts_est -qlow*ts_se),na.rm=T)
    ci_length_qt = apply((qhigh-qlow)*ts_se,1,median,na.rm=T)/(2*crit_norm*sd(mle_est))
    
    
    #Unfeasible bias-correction for CIs
    adjust_mle = 1/(ci_length_mle_est)
    adjust_est = 1/(ci_length_est/ci_length_true)  
    
    cov_mle_art =  mean((ff>=mle_est-adjust_mle*crit_norm*mle_se)&(ff<=mle_est+adjust_mle*crit_norm*mle_se))
    cov_lmoment_art = rowMeans((ff>=ts_est-crit_norm*ts_se*adjust_est)&(ff<=ts_est+crit_norm*ts_se*adjust_est),na.rm=T)
    
    #Calculate an (unfeasible) correction for the MLE that restores coverage to the nominal level
    qt_mle = quantile((mle_est-ff)/mle_se, c((1-coverage)/2, 1 - (1-coverage)/2))
    inflate = (diff(qt_mle)/(2*crit_norm))*(median(mle_se)/sd(mle_est))
  
    
    

    
    gridL = length(true.par) + 1:length(cov_lmoment) - 1
    
    
    if(!dir.exists(paste("results/",mc.name,"/coverage/",sep="")))
      dir.create(paste("results/",mc.name,"/coverage/",sep=""),recursive = T)
    
    
    
    #Coverage results
    pdf(paste((paste("results/",mc.name,"/coverage/coverage_",mc.name, "_N",ss,"_tau", tau, "_true.pdf",sep=""))),width= 6, height=6)
    
    plot(gridL, cov_lmoment_true, lty = 1,type = "l", col = 'blue',ylab = 'Coverage', xlab = "R",
         main = bquote(tau~"="~.(tau)), ylim = c(0.5,1), xlim=c(min(gridL),max(gridL)), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, lwd = 2)
    
    lines(gridL, rep(cov_mle_true,length(gridL)), col= 'red')
    abline(h=0.95, col = 'black')
    grid(col='gray')
    legend("bottomright", c("MLE", "Càglàd TS"),
           col = c("red", "blue"), lty = c(1,1))
    

    dev.off()
    
    
    pdf(paste((paste("results/",mc.name,"/coverage/coverage_",mc.name, "_N",ss,"_tau", tau, "_est.pdf",sep=""))),width= 6, height=6)
    
    plot(gridL, cov_lmoment, lty = 1,type = "l", col = 'blue',ylab = 'Coverage', xlab = "R",
         main = bquote(tau~"="~.(tau)), ylim = c(0.5,1), xlim=c(min(gridL),max(gridL)), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, lwd = 2)
    
    lines(gridL, rep(cov_mle,length(gridL)), col= 'red')
    abline(h=0.95, col = 'black')
    grid(col='gray')
    legend("bottomright", c("MLE", "Càglàd TS"),
           col = c("red", "blue"), lty = c(1,1))
    
    
    dev.off()
    
    
    pdf(paste((paste("results/",mc.name,"/coverage/length_",mc.name, "_N",ss,"_tau", tau, "_true.pdf",sep=""))),width= 6, height=6)
    
    plot(gridL, ci_length_true, lty = 1,type = "l", col = 'blue',ylab = 'Relative Average Length', xlab = "R",
         main = bquote(tau~"="~.(tau)), ylim = c(round(min(c(ci_length_true,0.5)),digits=2)-0.1,round(max(c(ci_length_true,1)),digits=2)+0.1), xlim=c(min(gridL),max(gridL)), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, lwd = 2)
    
    
    lines(gridL, rep(1,length(gridL)), col = 'red')
    
    grid(col='gray')
    legend("bottomright", c("MLE", "Càglàd TS"),
           col = c("red", "blue"), lty = c(1,1))
    
    
    dev.off()
    
    
    pdf(paste((paste("results/",mc.name,"/coverage/length_",mc.name, "_N",ss,"_tau", tau, "_est.pdf",sep=""))),width= 6, height=6)
    
    plot(gridL, ci_length_est, lty = 1,type = "l", col = 'blue',ylab = 'Relative Average length', xlab = "R",
         main = bquote(tau~"="~.(tau)), ylim =c(round(min(c(ci_length_est,0.5)),digits=2)-0.1,round(max(c(ci_length_est,1)),digits=2)+0.1), xlim=c(min(gridL),max(gridL)), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, lwd = 2)
    
    
    lines(gridL, rep(ci_length_mle_est,length(gridL)), col = 'red')
    grid(col='gray')
    legend("bottomright", c("MLE", "Càglàd TS"),
           col = c("red", "blue"), lty = c(1,1))

    dev.off()
    
    
    

    
    if((tau>0.5)&(ss<500))
    {
      
      #Coverage results
      pdf(paste((paste("results/",mc.name,"/coverage/coverage_",mc.name, "_N",ss,"_tau", tau, "_art_corr.pdf",sep=""))),width= 6, height=6)
      
      plot(gridL, cov_lmoment_art, lty = 1,type = "l", col = 'blue',ylab = 'Coverage', xlab = "R",
           main = bquote(tau~"="~.(tau)), ylim = c(0.5,1), xlim=c(min(gridL),max(gridL)), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, lwd = 2)
      
      lines(gridL, rep(cov_mle_art,length(gridL)), col= 'red')
      abline(h=0.95, col = 'black')
      grid(col='gray')
      legend("bottomright", c("MLE (rescaled)", "Càglàd TS (rescaled)"),
             col = c("red", "blue"), lty = c(1,1))
      
      
      dev.off()
      
      #Coverage correction
      pdf(paste((paste("results/",mc.name,"/coverage/coverage_",mc.name, "_N",ss,"_tau", tau, "_correct.pdf",sep=""))),width= 6, height=6)
      
      plot(gridL, cov_lmoment_qt, lty = 1,type = "l", col = 'blue',ylab = 'Coverage', xlab = "R",
           main = bquote(tau~"="~.(tau)), ylim = c(0.5,1), xlim=c(min(gridL),max(gridL)), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, lwd = 2)
      
      lines(gridL, rep(cov_mle,length(gridL)), col= 'red')
      abline(h=0.95, col = 'black')
      
      grid(col='gray')
      legend("bottomright", c("MLE (uncorrected)", "Càglàd TS (corrected)"),
             col = c("red", "blue"), lty = c(1,1))
      
      
      dev.off()
      
      
      #Coverage correction
      pdf(paste((paste("results/",mc.name,"/coverage/length_",mc.name, "_N",ss,"_tau", tau, "_correct.pdf",sep=""))),width= 6, height=6)
      
      plot(gridL, ci_length_qt, lty = 1,type = "l", col = 'blue',ylab = 'Relative Average length', xlab = "R",
           main = bquote(tau~"="~.(tau)), ylim =c(round(min(c(ci_length_qt,0.5,inflate,ci_length_mle_est)),digits=2)-0.1,round(max(c(ci_length_qt,1, inflate, ci_length_mle_est)),digits=2)+0.1), xlim=c(min(gridL),max(gridL)), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, lwd = 2)
      
      
      lines(gridL, rep(ci_length_mle_est,length(gridL)), col = 'red')
      lines(gridL, rep(inflate, length(gridL)), col = 'red', lty=2)
      
      grid(col='gray')
      legend("bottomright", c("MLE (uncorrected)","MLE (unfeasibly corrected)", "Càglàd TS (corrected)"),
             col = c("red","red", "blue"), lty = c(1,2,1))
      
      dev.off()
      
    }
    
   
  }
  
}
