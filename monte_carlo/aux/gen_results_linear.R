library(xtable)
source("../methods/estimation.R")

table = c()

table_tilt = c()
table_trim = c()

#Matrix power 
"%^%" <- function(S, power) 
  with(eigen(S), vectors %*% (diag(values^power)) %*% t(vectors)) 
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
    
    saveRDS(keeper, paste(mc.name,"_mc/keeper_",mc.name, "_N",ss, ".RDS",sep=""))
    
  }else print("Everything ok!")
  
  
  
  mle = Reduce("+",lapply(modelo, function(x){ zz = x$mle$par -true.par
               
               zz%*%t(zz)}))/length(modelo)
  
  lmoments.fs.cag.mat = sapply(modelo, function(x){
    lapply(x$caglad_fs, function(y){
      zz = y$fs$par-true.par
    zz%*%t(zz)})})
  
  lmoments.fs.cag = lapply(1:nrow(lmoments.fs.cag.mat), function(x) Reduce('+',lmoments.fs.cag.mat[x,])/length(modelo))
  
  
  lmoments.ss.cag.mat = sapply(modelo, function(x){
    lapply(x$caglad_ss, function(y){
      zz = y$ss$par-true.par
      zz%*%t(zz)})})
  
  lmoments.ss.cag = lapply(1:nrow(lmoments.ss.cag.mat), function(x) Reduce('+',lmoments.ss.cag.mat[x,])/length(modelo))
  
  #Unbiased
  lmoments.fs.un.mat = sapply(modelo, function(x){
    lapply(x$unbiased_fs, function(y){
      zz = y$fs$par-true.par
      zz%*%t(zz)})})
  
  lmoments.fs.un = lapply(1:nrow(lmoments.fs.un.mat), function(x) Reduce('+',lmoments.fs.un.mat[x,])/length(modelo))
  
  #unbiased
  lmoments.ss.un.mat = sapply(modelo, function(x){
    lapply(x$unbiased_ss, function(y){
      zz = y$ss$par-true.par
      zz%*%t(zz)})})
  
  lmoments.ss.un = lapply(1:nrow(lmoments.ss.un.mat), function(x) Reduce('+',lmoments.ss.un.mat[x,])/length(modelo))


  
  B_mat = mle%^%(-1/2)
  
  
  for(cases in c("lmoments.fs.cag","lmoments.ss.cag", "lmoments.fs.un", "lmoments.ss.un"))
    assign(cases, sapply(get(cases), function(x){
      
      egg = eigen(B_mat%*%x%*%B_mat)
      
      sqrt(egg$values[rev(c(1,length(true.par)))])
    }))
  
  
    
    
    line.best = format(round(c(min(lmoments.fs.cag[1,]),min(lmoments.ss.cag[1,]), min(lmoments.fs.un[1,]), min(lmoments.ss.un[1,])),digits=3),scientific=F)
    line.which.best = c(which.min(lmoments.fs.cag[1,]),which.min(lmoments.ss.cag[1,]), which.min(lmoments.fs.un[1,]), which.min(lmoments.ss.un[1,])) + length(true.par) - 1
    
    line.worst = format(round(c(min(lmoments.fs.cag[2,]),min(lmoments.ss.cag[2,]), min(lmoments.fs.un[2,]), min(lmoments.ss.un[2,])),digits=3),scientific=F)
    line.which.worst = c(which.min(lmoments.fs.cag[2,]),which.min(lmoments.ss.cag[2,]), which.min(lmoments.fs.un[2,]), which.min(lmoments.ss.un[2,])) + length(true.par) - 1
    
    line.which.best = paste(" (",line.which.best, ")",sep="")
    line.which.worst = paste(" (",line.which.worst, ")",sep="")
    
    table = cbind(table,
                  cbind(as.vector(rbind(line.best,line.which.best)),as.vector(rbind(line.worst,line.which.worst))))
    

  
}
#Rearranging so most favourable first and least favourable later
table = cbind(table[,seq(1,length.out=length(SampleSize),by=2)],table[,seq(2,length.out=length(SampleSize),by=2)])

table = cbind(as.vector(rbind(c("Càglàd FS", "Càglàd TS", "Unbiased FS", "Unbiased TS"),"")),table)

titulo = c("",rep(paste("T=", SampleSize), 2))

table = rbind(titulo,table)

alignment = c("|l|","|l|")
alignment_base = c(rep("c", length(SampleSize)-1),"c|")

alignment = c(alignment, rep(alignment_base,2))

header =paste("&", paste(paste("\\multicolumn{",length(SampleSize),"}","{c|}{"),c("Most favourable $\\delta$","Least favourable $\\delta$"),"}",
                         collapse = "&"), "\\\\")

print(xtable(table,caption = paste(toupper(mc.name),": maximal and minimal relative RMSE for linear combintations of parameters (MSE-minimising choice of $R$)"),
             label = paste(mc.name, "_table_linear",sep=""), align = alignment), type = "latex",
      file = paste("results/",mc.name, "_table_linear.tex",sep=""), table.placement = "H",
      caption.placement = "top", include.rownames = F, include.colnames = F,
      sanitize.text.function = identity, 
      add.to.row = list("pos"=list(0),
                        "command" = header),
      hline.after = c(-1,1))



