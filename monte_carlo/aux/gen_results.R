library(xtable)
source("../methods/estimation.R")

table = c()

for(ss in SampleSize)
{

modelo = readRDS(paste(mc.name, "_N",ss, ".RDS",sep=""))


#In this step, we rerun estimators which failed (generated errors) due to convergence issues

#Loads data
yMat = readRDS( paste("sample_",mc.name,"_N",ss, ".RDS",sep=""))


#Checks for caglad
val = which(do.call(rbind,lapply(modelo, function(x){sapply(x$caglad_ss, function(y){length(y)!=2})})),arr.ind=T)
if(nrow(val)>0)
  for(j in 1:nrow(val))
    modelo[[val[j,1]]]$caglad_ss[[val[j,2]]] = lmoment.est(yMat[,val[j,1]], modelo[[val[j,1]]]$caglad_fs[[1]]$fs$par,val[j,2] + length(true.par)-1 , lmoment.analytic = lmoment.analytic, quantile.func = quantile.function,
                                                           density = density.function, grid.length = 2000, par.first.step = modelo[[val[j,1]]]$caglad_fs[[1]]$fs$par, control = list( "maxit"=500))

#Checks for unbiased
val = which(do.call(rbind,lapply(modelo, function(x){sapply(x$unbiased_ss, function(y){length(y)!=2})})),arr.ind=T)
if(nrow(val)>0)
  for(j in 1:nrow(val))
    modelo[[val[j,1]]]$unbiased_ss[[val[j,2]]] = lmoment.est(yMat[,val[j,1]], modelo[[val[j,1]]]$unbiased_fs[[1]]$fs$par,val[j,2] + length(true.par)-1 , lmoment.analytic = lmoment.analytic, quantile.func = quantile.function,
                                                             lmoment.est = "unbiased", density = density.function, grid.length = 2000, par.first.step = modelo[[val[j,1]]]$unbiased_fs[[1]]$fs$par, control = list( "maxit"=500))


#For the GEV MC, recalculate FS estimator that diverged
if((mc.name == "gev")&(ss==100))
{
  print(modelo[[907]]$caglad_fs[[2]])
  
  modelo[[907]]$caglad_fs[[2]] = lmoment.est(yMat[,907], modelo[[907]]$caglad_fs[[1]]$fs$par, 2+length(true.par)-1, lmoment.analytic = lmoment.analytic, quantile.func = quantile.function,
                                           density = density.function, grid.length = 2000, weight.matrix = "first", control = list( "maxit"=500))

}

for(tau in tau.seq)
{
  ff = quantile.function(tau, true.par)
  
  
  
  mle = mean(sapply(modelo, function(x){mean((quantile.function(tau,x$mle$par)-ff)^2)}))
  
  lmoments.fs.cag = colMeans(do.call(rbind,lapply(modelo, function(x){sapply(x$caglad_fs, function(y){mean((quantile.function(tau,y$fs$par)-ff)^2)})})))

  lmoments.ss.cag = colMeans(do.call(rbind,lapply(modelo, function(x){sapply(x$caglad_ss, function(y){mean((quantile.function(tau,y$ss$par)-ff)^2)})})))
  
  lmoments.fs.un = colMeans(do.call(rbind,lapply(modelo, function(x){sapply(x$unbiased_fs, function(y){mean((quantile.function(tau,y$fs$par)-ff)^2)})})))
  
  lmoments.ss.un = colMeans(do.call(rbind,lapply(modelo, function(x){sapply(x$unbiased_ss, function(y){mean((quantile.function(tau,y$ss$par)-ff)^2)})})))
  
  
  lmoments.fs.cag = sqrt(lmoments.fs.cag)/sqrt(mle) 
  lmoments.ss.cag = sqrt(lmoments.ss.cag)/sqrt(mle) 
  lmoments.fs.un = sqrt(lmoments.fs.un)/sqrt(mle) 
  lmoments.ss.un = sqrt(lmoments.ss.un)/sqrt(mle) 
  
  
  gridY = c(min(min(lmoments.ss.cag,lmoments.fs.cag, lmoments.fs.un, lmoments.ss.un),1), max(max(lmoments.ss.cag,lmoments.fs.cag, lmoments.fs.un, lmoments.ss.un),1))
  gridL = length(true.par) + 1:length(lmoments.fs.cag) - 1
  
  
  if(!dir.exists(paste("results/",mc.name,sep="")))
    dir.create(paste("results/",mc.name,sep=""),recursive = T)
  
  pdf(paste((paste("results/",mc.name,"/",mc.name, "_N",ss,"_tau", tau, ".pdf",sep=""))),width= 6, height=6)
  
  plot(gridL, lmoments.fs.cag, lty = 1,type = "l", ylim = gridY, col = 'blue',ylab = 'Relative RMSE', xlab = "R",
       main = bquote(tau~"="~.(tau)), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, lwd = 2)
  lines(gridL, lmoments.ss.cag, lty = 2, col='blue', lwd = 2)
  
  lines(gridL, lmoments.fs.un, lty = 1, col = 'red', lwd = 2 )
  lines(gridL, lmoments.ss.un, lty = 2, col = 'red', lwd = 2)
  
  legend("topleft", c("Càglàd FS", "Càglàd TS", "Unbiased FS", "Unbiased TS"),
         col = c("blue", "blue", "red", "red"), lty = c(1,2,1,2), lwd = c(2,2,2,2), bty = "n", cex = 1.5)
  
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