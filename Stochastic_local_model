rm(list = ls())
library(deSolve)

Logistic_mod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    if(N<0){
      N <- 0
    }
    dy <- (r + rnorm(1,0,sigma))*N*(1-N/K)
    
    return(list(dy))
  })
}

rep_time <- 0
survive <- 0
extinct <- 0
repeat{
  rep_time <- rep_time + 1
  pars <- c(r = 1, K = 1, sigma = 7)
  times <- seq(0, 50, by = 0.1)
  ini <- c(N = 0.5)
  out <- ode(ini, times, Logistic_mod, pars, method = "ode45", atol = 0.01, rtol = 0.01)
  #plot(out)
  if(min(out[,2], na.rm = TRUE) >= 0){
    survive <- survive + 1
  } else{
    extinct <- extinct + 1
  }
  if(rep_time>=100){
    break
  }
}

rep_time <- 0
coexist <- 0
exclusion <- 0
collapse <- 0

GNMmod <- function(Time, n, Pars) {
  with( as.list(Pars), {
    n[n < 0] <- 0
    dn <- (r + rnorm(3,0,sigma))*n*(1 - n + bij%*%n - aij%*%n^2)
    return(list(dn))
  })
}
ini <- c(0.5, 0.5, 0.5)
aij <- matrix(c(0,2,0,0,0,3,3,0,0),ncol = 3)
bij <- matrix(c(0,2,0,0,0,3,3,0,0),ncol = 3)
r <- c(1,1)
pars <- c(r = r, aij = aij, bij = bij, sigma = 0.8)
times <- seq(0, 200, by = 0.1)
out <- ode(ini, times, GNMmod, pars, method = "ode45", atol = 0.01, rtol = 0.01)
plot(out)


