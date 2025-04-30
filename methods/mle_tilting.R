empirically.tilted <- function(y, f_density, par0, delta, method = 'BFGS', control = list('maxit'=500))
{
  objective <- function(par) {
    dens = f_density(y, par)
    
    root.find <- function(x) {
      x*sum(dens^x*log(dens))/sum(dens^x) - log(sum(dens^x)) + log(length(y)) - delta
    }

    c2 <- uniroot(root.find, lower  = 0, upper =100)$root
    
    c1 = 1/(sum(dens^c2))
    return(-sum(c1*dens^c2*log(dens)))
  }
  
  optim(par0, objective, method = method, control = control)
}