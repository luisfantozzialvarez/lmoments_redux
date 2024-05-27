rm(list=ls())

setwd("~/Documents/GitHub/lmoments_redux/application/")

julia_path = '/Applications/Julia-1.9.app/Contents/Resources/julia/bin'

library(plm)
library(grf)

dados = read.csv("dataset_ridesharing.csv")



semanas = unique(dados$week)


early = dados[dados$week!="2018-09-24",]
late = dados[dados$week=="2018-09-24",]


set.seed(123)

source("specifications_application.R")
source("../methods/estimation.R")
source("../methods/selection.R")
source("../methods/semiparametric_model.R")


#Step 1: preliminary Pareto estimation for raw data

fs = lmoment.est(late$expenses, c(1,1), 2, lmoment.analytic = lmoment.analytic.gpd, quantile.func = quantile.function.gpd,
                 density = density.function.gpd, grid.length = 2000, weight.matrix = "first", control = list( "maxit"=500))

tau.seq = seq(0.1,0.9, by = c(0.1))

tabua = lmoment.select(late$expenses, fs$fs$par, 150, orthogonal = F, uvalues = tau.seq, lmoment.analytic = lmoment.analytic.gpd, quantile.func=quantile.function.gpd, density.function = density.function.gpd,
                       lmoment.est = "caglad", grid.length = 2000, Nsim = 1000,
                       control = list( "maxit"=500),mc.cores= detectCores() )

L= median(tabua$Lvals)

ts = lmoment.est(late$expenses, fs$fs$par,L, lmoment.analytic = lmoment.analytic.gpd, quantile.func = quantile.function.gpd,
                 density = density.function.gpd, grid.length = 2000, weight.matrix = "par", control = list( "maxit"=500),
                 vcov=T)


#Step 2: estimating the semiparametric mixture model


#Running preliminary estimator
y_pos = which(colnames(dados)=="dex")
x_pos = 4:484
w_pos = which(colnames(dados)=="ldex1")
s_pos = which(colnames(dados)=="lex2")

#First step, just to estimate residuals
model <- lmoment.semiparametric(early, late, y_pos, x_pos, w_pos, s_pos, c(1,1,1), L, grid_inv = 1000, 
                                lmoment.analytic = lmoment.analytic.gev, quantile.func = quantile.function.gev,
                                density.function = density.function.gev,   control = list( "maxit"=500))




#Now, finding mixture. We run estimator from different starting mixing points
model.twoside = model$model

grid.test <- lapply(seq(0.01,0.99,0.01), function(u) {
  
  print(u)
  
  model.lower = model$model
  model.lower$res_oos = model.lower$res_oos[order(model.lower$res_oos)]
  model.lower$corr = model.lower$corr[order(model.lower$corr)]
  
  
  model.lower$res_oos = model.lower$res_oos[1:ceil(u*length(model.lower$res_oos))]
  model.lower$corr = model.lower$corr[1:ceil(u*length(model.lower$corr))]
  
  
  model.upper = model$model
  model.upper$res_oos = model.upper$res_oos[order(model.upper$res_oos)]
  model.upper$corr = model.upper$corr[order(model.upper$corr)]
  
  
  model.upper$res_oos = model.upper$res_oos[(ceil(u*length(model.upper$res_oos))+1):(length(model.upper$res_oos))]
  model.upper$corr = model.upper$corr[(ceil(u*length(model.upper$corr))+1):(length(model.upper$corr))]
  

  lower <- lapply(c(-1,1), function(x){
    lmoment.semiparametric(early, late, y_pos, x_pos, w_pos, s_pos,  model$est$par, L, grid_inv = 1000, 
                           lmoment.analytic = lmoment.analytic.gev, quantile.func = quantile.function.gev,
                           density.function = density.function.gev,  
                           model = list('res_oos' = x*model.lower$res_oos,
                                        'corr' = x*model.lower$corr), control = list( "maxit"=1000)
    )
  })
  
  upper <- lapply(c(-1,1), function(x){
    lmoment.semiparametric(early, late, y_pos, x_pos, w_pos, s_pos,  model$est$par, L, grid_inv = 1000, 
                           lmoment.analytic = lmoment.analytic.gev, quantile.func = quantile.function.gev,
                           density.function = density.function.gev,   
                           model = list('res_oos' = x*model.upper$res_oos,
                                        'corr' = x*model.upper$corr), control =list( "maxit"=1000)
    )
  })
  
  
  par_guess.1 =  c(lower[[which.min(sapply(lower, function(x) x$est$value))]]$est$par,u,upper[[which.min(sapply(upper, function(x) x$est$value))]]$est$par)
  
  fs.1 <- tryCatch({lmoment.semiparametric(early, late, y_pos, x_pos, w_pos, s_pos, par_guess.1,
                                         L, grid_inv = 1000, 
                                         lmoment.analytic = lmoment.analytic.gev.mix, quantile.func =  quantile.function.gev.mix,
                                         density.function = density.function.gev.mix, bounds = list('lower'=c(-Inf,0,-Inf,0,-Inf,0,-Inf),'upper'=c(Inf,Inf,Inf,1,Inf,Inf,Inf)),
                                         model = model.twoside, control = list( "maxit"=1000))},error = function(e){ print(e); return(NULL)})
                                         
  model.tail.1 <- tryCatch({lmoment.semiparametric(early, late, y_pos, x_pos, w_pos, s_pos, par_guess.1,
                                                 L, grid_inv = 1000, 
                                                 lmoment.analytic = lmoment.analytic.gev.mix, quantile.func =  quantile.function.gev.mix,
                                                 density.function = density.function.gev.mix, 
                                                 model = model.twoside, control = list( "maxit"=1000), first.step = fs.1$est)}, error = function(e) NULL)
  
  par_guess.2 =  c(upper[[which.min(sapply(upper, function(x) x$est$value))]]$est$par, 1-u, lower[[which.min(sapply(lower, function(x) x$est$value))]]$est$par)
  
  fs.2 <- tryCatch({lmoment.semiparametric(early, late, y_pos, x_pos, w_pos, s_pos, par_guess.2,
                                           L, grid_inv = 1000, 
                                           lmoment.analytic = lmoment.analytic.gev.mix, quantile.func =  quantile.function.gev.mix,
                                           density.function = density.function.gev.mix, bounds = list('lower'=c(-Inf,0,-Inf,0,-Inf,0,-Inf),'upper'=c(Inf,Inf,Inf,1,Inf,Inf,Inf)),
                                           model = model.twoside, control = list( "maxit"=1000))},error = function(e){ print(e); return(NULL)})
  
  model.tail.2 <- tryCatch({lmoment.semiparametric(early, late, y_pos, x_pos, w_pos, s_pos, par_guess.2,
                                                   L, grid_inv = 1000, 
                                                   lmoment.analytic = lmoment.analytic.gev.mix, quantile.func =  quantile.function.gev.mix,
                                                   density.function = density.function.gev.mix, 
                                                   model = model.twoside, control = list( "maxit"=1000), first.step = fs.2$est)}, error = function(e) NULL)
  
  
   return(list('tail1' = model.tail.1, 'tail2' = model.tail.2, 'par.1' = par_guess.1, 'par.2' = par_guess.2))
})

save.image('intermmediate.RDS')
position = sapply(grid.test, function(x) c(ifelse(is.null(x[[1]]), Inf, x[[1]]$est$value), ifelse(is.null(x[[2]]), Inf, x[[2]]$est$value)))
line_min = which.min(apply(position,2 ,min))

col_min = which.min(position[,line_min])


#This will be the mixture model
model.tail = grid.test[[line_min]][[col_min]]

prev.tail <-model.tail
#We will now apply a continuously-updating procedure to improve the weighting matrix estimator
list.tail <- list(model.tail)
maxit=30
for(j in 1:20)
  if(!is.null(list.tail[[j]]))
    list.tail[[j+1]] <- tryCatch({lmoment.semiparametric(early, late, y_pos, x_pos, w_pos, s_pos, par_guess.2,
                                                 L, grid_inv = 1000, 
                                                 lmoment.analytic = lmoment.analytic.gev.mix, quantile.func =  quantile.function.gev.mix,
                                                 density.function = density.function.gev.mix, 
                                                 model = model.twoside, control = list( "maxit"=1000), first.step = list.tail[[j]]$est)}, error = function(e) NULL)


model.tail= list.tail[[1+which.min(sapply(list.tail[2:(maxit+1)], function(x) ifelse(is.null(x), Inf, x$est$value)))]]


#3. Generating results
grid_x = seq(min(late$expenses),max(late$expenses),by=0.1)
dens <- function(par) density.function.gpd(grid_x, par)
jacob_gpd = ad_jacobian(dens,ts$ss$par)

grid_e = seq(min(model$model$res_oos),to=max(model$model$res_oos),0.1)
dens_e <- function(par) density.function.gev.mix(grid_e,par)
jacob_e = jacob.density.function.gev.mix(grid_e,model.tail$est$par)


#Results
#Setting seed again
set.seed(123)

#3.1 Raw data

#Overidentifying restrictions test
pvalue = 1 - pchisq(ts$ss$value*nrow(late), L-length(ts$ss$par))


Nsims = 5000

vcov.gpd = ts$vcov/ts$N

sim_crit = jacob_gpd%*%t(chol(vcov.gpd))%*%matrix(rnorm(length(ts$ss$par)*Nsims),length(ts$ss$par))

diag_center = sqrt(diag(jacob_gpd%*%vcov.gpd%*%t(jacob_gpd)))
sim_crit = sim_crit/diag_center

vlh = apply(sim_crit,2, function(x)max(abs(x)))

lower_upper = cbind(dens(ts$ss$par) - diag_center*quantile(vlh,0.95), dens(ts$ss$par) + diag_center*quantile(vlh,0.95))


curve.f <- function(x) density.function.gpd(x, ts$ss$par)


pdf('expenses.pdf', height = 5, width=5)
hist(late$expenses,freq=F,breaks = seq(min(10*ceiling(late$expenses)%/%10),max(11*floor(late$expenses)%/%10),10), xlab = 'Expenditures (BRL)',
     main = 'Weekly expenditure')
polygon(c(grid_x, rev(grid_x)), c(lower_upper[,1], rev(lower_upper[,2])), col = "#0000FF50",border=NA)

curve(curve.f, from = min(late$expenses), to = max(late$expenses), col = 'blue', add=T)
dev.off()



#3.2 Model errors

#Overidentifying restrictions test
pvalue = 1 - pchisq(model.tail$est$value*model.tail$N, L-length(model.tail$est$par))


vcov.gev = model.tail$vcov/model.tail$N

sim_crit = jacob_e%*%t(chol(vcov.gev))%*%matrix(rnorm(length(model.tail$est$par)*Nsims),nrow=length(model.tail$est$par))

diag_center = sqrt(diag(jacob_e%*%vcov.gev%*%t(jacob_e)))
sim_crit = sim_crit/diag_center

vlh = apply(sim_crit,2, function(x)max(abs(x)))

lower_upper = cbind(dens_e(model.tail$est$par) - diag_center*quantile(vlh,0.95), dens_e(model.tail$est$par) + diag_center*quantile(vlh,0.95))


curve.f <- function(x) density.function.gev.mix(x, model.tail$est$par)

pdf('res.pdf', height = 5, width=5)

hist(model$model$res_oos,freq=F,breaks = seq(min(10*ceiling(model$model$res_oos)%/%10),max(11*floor(model$model$res_oos)%/%10),10), xlab = 'Expenditures (BRL)',
     main =  expression(Delta~hat(e)[it]))

polygon(c(grid_e, rev(grid_e)), c(lower_upper[,1], rev(lower_upper[,2])), col = "#0000FF50",border=NA)

curve(curve.f,min(model$model$res_oos),max(model$model$res_oos),col = 'blue', add=T)
dev.off()

#99 Percent CI

write.csv(cbind(model.tail$est$par, model.tail$est$par - qnorm(1-0.01/2)*sqrt(diag(vcov.gev)),model.tail$est$par +  qnorm(1-0.01/2)*sqrt(diag(vcov.gev))),'coefs.csv')

save.image('final.RData')



