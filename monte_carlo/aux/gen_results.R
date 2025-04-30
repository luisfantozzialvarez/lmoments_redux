library(xtable)
source("../methods/estimation.R")

table = c()

table_tilt = c()
table_trim = c()
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
  
  #Saving nonconvergent draws for later use in other tables
  saveRDS(keeper, paste(mc.name,"_mc/keeper_",mc.name, "_N",ss, ".RDS",sep=""))
  
  }else print("Everything ok!")





for(tau in tau.seq)
{
  ff = quantile.function(tau, true.par)
  
  
  
  mle = mean(sapply(modelo, function(x){mean((quantile.function(tau,x$mle$par)-ff)^2)}))
  
  lmoments.fs.cag = colMeans(do.call(rbind,lapply(modelo, function(x){sapply(x$caglad_fs, function(y){mean((quantile.function(tau,y$fs$par)-ff)^2)})})))

  lmoments.ss.cag = colMeans(do.call(rbind,lapply(modelo, function(x){sapply(x$caglad_ss, function(y){mean((quantile.function(tau,y$ss$par)-ff)^2)})})))
  
  lmoments.fs.un = colMeans(do.call(rbind,lapply(modelo, function(x){sapply(x$unbiased_fs, function(y){mean((quantile.function(tau,y$fs$par)-ff)^2)})})))
  
  lmoments.ss.un = colMeans(do.call(rbind,lapply(modelo, function(x){sapply(x$unbiased_ss, function(y){mean((quantile.function(tau,y$ss$par)-ff)^2)})})))
  
  trimmed_mle = colMeans(do.call(rbind,lapply(modelo, function(x){sapply(x$trimmed_mle, function(y){if(!is.null(y$par)) mean((quantile.function(tau,y$par)-ff)^2) else NA})})),na.rm=T)
  
  tilted_mle = colMeans(do.call(rbind,lapply(modelo, function(x){sapply(x$tilted_mle, function(y){if(!is.null(y$par)) mean((quantile.function(tau,y$par)-ff)^2) else NA})})),na.rm=T)
  
  trimmed_mle_prop = colMeans(is.na(do.call(rbind,lapply(modelo, function(x){sapply(x$trimmed_mle, function(y){if(!is.null(y$par)) mean((quantile.function(tau,y$par)-ff)^2) else NA})}))))
  tilted_mle_prop = colMeans(is.na(do.call(rbind,lapply(modelo, function(x){sapply(x$tilted_mle, function(y){if(!is.null(y$par)) mean((quantile.function(tau,y$par)-ff)^2) else NA})}))))
  
  
  lmoments.fs.cag = sqrt(lmoments.fs.cag)/sqrt(mle) 
  lmoments.ss.cag = sqrt(lmoments.ss.cag)/sqrt(mle) 
  lmoments.fs.un = sqrt(lmoments.fs.un)/sqrt(mle) 
  lmoments.ss.un = sqrt(lmoments.ss.un)/sqrt(mle) 
  
  trimmed_mle = sqrt(trimmed_mle/mle)
  tilted_mle = sqrt(tilted_mle/mle)
  
  gridY = c(min(min(lmoments.ss.cag,lmoments.fs.cag, lmoments.fs.un, lmoments.ss.un),1), max(max(lmoments.ss.cag,lmoments.fs.cag, lmoments.fs.un, lmoments.ss.un),1))
  gridL = length(true.par) + 1:length(lmoments.fs.cag) - 1
  
  
  if(!dir.exists(paste("results/",mc.name,sep="")))
    dir.create(paste("results/",mc.name,sep=""),recursive = T)
  
  pdf(paste((paste("results/",mc.name,"/",mc.name, "_N",ss,"_tau", tau, ".pdf",sep=""))),width= 6, height=6)
  
  plot(gridL, lmoments.fs.cag, lty = 1,type = "l", ylim = gridY, col = 'blue',ylab = 'Relative RMSE', xlab = "R",
       main = bquote(tau~"="~.(tau)), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, lwd = 2)
  lines(gridL, lmoments.ss.cag, lty = 6, col='blue', lwd = 2)
  
  lines(gridL, lmoments.fs.un, lty = 5, col = 'red', lwd = 2 )
  lines(gridL, lmoments.ss.un, lty = 3, col = 'red', lwd = 3)
  
  grid(col='gray')
  
  legend("topleft", c("Càglàd FS", "Càglàd TS", "Unbiased FS", "Unbiased TS"),
         col = c("blue", "blue", "red", "red"), lty = c(1,6,5,3), lwd = c(2,2,2,3), bty = "n", cex = 1.5)
  
  abline(h = 1,col='black')
  dev.off()
  
  if(!dir.exists(paste(paper_path,"/plots/",mc.name,sep="")))
    dir.create(paste(paper_path,"/plots/",mc.name,sep=""),recursive = T)
  
  file.copy(paste((paste("results/",mc.name,"/",mc.name, "_N",ss,"_tau", tau, ".pdf",sep=""))),
            paste(paper_path,"/plots/",mc.name,"/",mc.name, "_N",ss,"_tau", tau, ".pdf",sep=""), overwrite = T )
  
  line = format(round(c(min(lmoments.fs.cag),min(lmoments.ss.cag), min(lmoments.fs.un), min(lmoments.ss.un)),digits=3),scientific=F)
  line.which = c(which.min(lmoments.fs.cag),which.min(lmoments.ss.cag), which.min(lmoments.fs.un), which.min(lmoments.ss.un)) + length(true.par) - 1
  
  line.which = paste(" (",line.which, ")",sep="")
  table = cbind(table,
    as.vector(rbind(line,line.which)))
  
  
  table_trim = cbind(table_trim, as.numeric(t(cbind(trimmed_mle,trimmed_mle_prop))))
  table_tilt = cbind(table_tilt, as.numeric(t(cbind(tilted_mle,tilted_mle_prop))))

}

}

table = cbind(as.vector(rbind(c("Càglàd FS", "Càglàd TS", "Unbiased FS", "Unbiased TS"),"")),table)

titulo = c("",rep(paste("$\\tau =",tau.seq,  "$"), length(SampleSize)))

table = rbind(titulo,table)

alignment = c("|l|","|l|")
alignment_base = c(rep("c", length(tau.seq)-1),"c|")
alignment = c(alignment, rep(alignment_base, length(SampleSize)))

header =paste("&", paste(paste("\\multicolumn{",length(tau.seq),"}","{c|}{$T ="),SampleSize,"$}",
                         collapse = "&"), "\\\\")

print(xtable(table,caption = paste(toupper(mc.name),": relative RMSE under MSE-minimising choice of $R$"),
       label = paste(mc.name, "_table_mle",sep=""), align = alignment), type = "latex",
      file = paste("results/",mc.name, "_table_mle.tex",sep=""), table.placement = "H",
      caption.placement = "top", include.rownames = F, include.colnames = F,
      sanitize.text.function = identity, scalebox = 0.6,
      add.to.row = list("pos"=list(0),
                        "command" = header),
      hline.after = c(-1,1, 1+seq(2,length.out = 4, by=2)))


if(!dir.exists(paste(paper_path,"/tables/",mc.name,sep="")))
  dir.create(paste(paper_path,"/tables/",mc.name,sep=""),recursive = T)

file.copy(paste("results/",mc.name, "_table_mle.tex",sep=""),
          paste(paper_path,"/tables/",mc.name,"/",mc.name, "_table_mle.tex",sep=""), overwrite = T )

tabla_adicional = format(round(rbind(table_trim,table_tilt),digits=3), scientific = F)

tabla_adicional = cbind( as.vector(rbind(as.vector(sapply(c("Trimming", "Tilting"),
                                      function(x)paste(x, paste(100*(1-tau.seq[tau.seq>0.5]),"\\%",sep="")))),"")),tabla_adicional)
  
tabla_adicional[c(2,4,6),] = ""

tabla_adicional[c(8,10,12),-1]= paste("(",round(as.numeric(tabla_adicional[c(8,10,12),-1]),digits=2)*100,"\\%)",sep="")

tabla_adicional = rbind(table[c(1,4:5),],tabla_adicional)


print(xtable(tabla_adicional,caption = paste(toupper(mc.name),": comparison with Trimmed and Tilted MLE methods"),
             label = paste(mc.name, "_table_additional",sep=""), align = alignment), type = "latex",
      file = paste("results/",mc.name, "_additional.tex",sep=""), table.placement = "H",
      caption.placement = "top", include.rownames = F, include.colnames = F,
      sanitize.text.function = identity, scalebox = 0.6,
      add.to.row = list("pos"=list(0),
                        "command" = header),
      hline.after = c(-1,1, seq(3,length.out = 3, by=2*length(tau.seq[tau.seq>0.5]))))
