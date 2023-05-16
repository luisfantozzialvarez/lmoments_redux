library(xtable)

table = c()

for(ss in SampleSize)
{


  modelo = readRDS(paste("select",mc.name, "_N",ss, ".RDS",sep=""))

  ff = quantile.function(tau.seq, true.par)
  
  mle = sqrt(rowMeans(sapply(modelo, function(x){(quantile.function(tau.seq,x$mle$par)-ff)^2})))
  
  caglad.fs = sqrt(rowMeans(sapply(modelo, function(x){(quantile.function(tau.seq,x$select_caglad$first.step$par)-ff)^2})))
  caglad.rmse = sqrt(rowMeans(sapply(modelo, function(x){
    
    sapply(1:length(tau.seq), function(j){(quantile.function(tau.seq[j],x$caglad_ss[[j]]$ss$par)- ff[j])^2})
   })))
  
  caglad.lasso = sqrt(rowMeans(sapply(modelo, function(x){(quantile.function(tau.seq,x$lasso_select$ss$par)-ff)^2})))
  caglad.postlasso = sqrt(rowMeans(sapply(modelo, function(x){(quantile.function(tau.seq,x$lasso_select$ps$par)-ff)^2})))
  
  caglad.fs = caglad.fs/mle
  caglad.rmse = caglad.rmse/mle
  caglad.lasso= caglad.lasso/mle
  caglad.postlasso= caglad.postlasso/mle
  
  medL = apply(sapply(modelo, function(x) x$select_caglad$Lvals),1,mean)
  lassoL = mean(sapply(modelo, function(x) x$lasso_select$kept) )
  
  line = format(round(rbind(caglad.fs,caglad.rmse, caglad.lasso, caglad.postlasso),digits=3), scientific = F)
  
  pos_change = (rbind(caglad.fs,caglad.rmse, caglad.lasso, caglad.postlasso)>10)
  
  line[pos_change] = "$>10$"
  
  line.L = round(rbind(rep(length(true.par),ncol(line)), medL, rep(lassoL, ncol(line)), rep(lassoL, ncol(line))),digits=2)
  
  line.L = paste(" (",line.L, ")",sep="")
  
  valor = matrix(as.vector(rbind(as.vector(line),as.vector(line.L))), nrow = nrow(line)*2)
  
  
 table = cbind(table,
                valor)
}

table = cbind(as.vector(rbind(c("FS", "TS RMSE", "TS Lasso", "TS Post-Lasso"),"")),table)

titulo = c("",rep(paste("$\\tau =",tau.seq,  "$"), length(SampleSize)))

table = rbind(titulo,table)

alignment = c("|l|","|l|")
alignment_base = c(rep("c", length(tau.seq)-1),"c|")
alignment = c(alignment, rep(alignment_base, length(SampleSize)))

header =paste("&", paste(paste("\\multicolumn{",length(tau.seq),"}","{c|}{$T ="),SampleSize,"$}",
                         collapse = "&"), "\\\\")

print(xtable(table,caption = paste(toupper(mc.name),": relative RMSE under different selection procedures"),
             label = paste(mc.name, "_table_select",sep=""), align = alignment), type = "latex",
      file = paste("results/",mc.name, "_table_select.tex",sep=""), table.placement = "H",
      caption.placement = "top", include.rownames = F, include.colnames = F,
      sanitize.text.function = identity, scalebox = 0.6,
      add.to.row = list("pos"=list(0),
                        "command" = header),
      hline.after = c(-1,1, 1+seq(2,length.out = 4, by=2)))

print(xtable(table[-(6:7),],caption = paste(toupper(mc.name),": relative RMSE under different selection procedures"),
             label = paste(mc.name, "_table_select_summary",sep=""), align = alignment), type = "latex",
      file = paste("results/",mc.name, "_table_select_summary.tex",sep=""), table.placement = "H",
      caption.placement = "top", include.rownames = F, include.colnames = F,
      sanitize.text.function = identity, scalebox = 0.6,
      add.to.row = list("pos"=list(0),
                        "command" = header),
      hline.after = c(-1,1, 1+seq(2,length.out = 3, by=2)))

if(!dir.exists(paste(paper_path,"/tables/",mc.name,sep="")))
  dir.create(paste(paper_path,"/tables/",mc.name,sep=""),recursive = T)

file.copy(paste("results/",mc.name, "_table_select.tex",sep=""),
          paste(paper_path,"/tables/",mc.name,"/",mc.name, "_table_select.tex",sep=""), overwrite = T )

file.copy(paste("results/",mc.name, "_table_select_summary.tex",sep=""),
          paste(paper_path,"/tables/",mc.name,"/",mc.name, "_table_select_summary.tex",sep=""), overwrite = T )