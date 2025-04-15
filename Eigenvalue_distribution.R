rm(list=ls())
path = "F:/AcademyData/Computer/rcode/ecology_model/GNM3_laplasian/local_dynamics/"
dir.create(path)
setwd(path)
library(deSolve)
library(ggplot2)
library(patchwork)
library(tidyr)
library(paletteer)
library(xlsx)
library(openxlsx)
library(Rlab)
library(numDeriv)

GNMmod <- function(Time, n, Pars) {
  with( as.list(Pars), {
    n[n < 0] <- 0
    dn <- r*n*(1 - n + bij%*%n - aij%*%(n^2))
    return(list(dn))
  })
}

GLVmod <- function(Time, n, Pars) {
  with( as.list(Pars), {
    dn <- r*n*(1 - aij%*%(n))
    return(list(dn))
  })
}

model_type = "GLV"
number_of_polulation = 2
U <- c(2,2)
iterations = 2000
nn <- number_of_polulation^2
times <- seq(0, iterations, by = 0.01)
Eigval <- c()
repeat_time = 0
repeat{
  repeat_time = repeat_time + 1
  aij <- matrix(rgamma(nn,shape = 1,scale = U[1]),
                nrow = number_of_polulation,ncol = number_of_polulation)
  bij <- matrix(rgamma(nn,shape = 1,scale = U[2]),
                nrow = number_of_polulation,ncol = number_of_polulation)
  diag(aij) <- 0
  diag(bij) <- 0
  r = runif(number_of_polulation,0.5,1)
  pars <- c(aij,bij,r)
  NN <- runif(number_of_polulation,0,1)
  tryCatch({
    out <- ode(NN, times, GLVmod, pars)
    J <- jacobian(function(x){as.vector(GLVmod(Time = 1, n = x, Pars = pars)[[1]])}, 
                  out[nrow(out),2:ncol(out)], method = "simple")
    Eigval <<- c(Eigval, eigen(J)$values[which(Re(eigen(J)$values) == max(Re(eigen(J)$values)))])
  }, error = function(e){repeat_time <<- repeat_time - 1})
  if(repeat_time >= 5000){
    break
  }
}

mydf1 <- data.frame(x = Re(Eigval), y = Im(Eigval))
myp1 <- ggplot(mydf1, aes(x = x, y = y)) +
  geom_point(color = "red", size = 1) +
  labs(
    title = "n = 2",
    x = "Real",
    y = "Imaginary"
  ) +
  theme_minimal()

number_of_polulation = 3
Eigval <- c()
repeat_time = 0
repeat{
  repeat_time = repeat_time + 1
  aij <- matrix(rgamma(nn,shape = 1,scale = U[1]),
                nrow = number_of_polulation,ncol = number_of_polulation)
  bij <- matrix(rgamma(nn,shape = 1,scale = U[2]),
                nrow = number_of_polulation,ncol = number_of_polulation)
  diag(aij) <- 0
  diag(bij) <- 0
  r = runif(number_of_polulation,0.5,1)
  pars <- c(aij,bij,r)
  NN <- runif(number_of_polulation,0,1)
  tryCatch({
    out <- ode(NN, times, GLVmod, pars)
    J <- jacobian(function(x){as.vector(GLVmod(Time = 1, n = x, Pars = pars)[[1]])}, 
                  out[nrow(out),2:ncol(out)], method = "simple")
    Eigval <<- c(Eigval, eigen(J)$values[which(Re(eigen(J)$values) == max(Re(eigen(J)$values)))])
  }, error = function(e){repeat_time <<- repeat_time - 1})
  if(repeat_time >= 5000){
    break
  }
}

mydf2 <- data.frame(x = Re(Eigval), y = Im(Eigval))
myp2 <- ggplot(mydf2, aes(x = x, y = y)) +
  geom_point(color = "red", size = 1) +
  labs(
    title = "n = 3",
    x = "Real",
    y = "Imaginary"
  ) +
  theme_minimal()

number_of_polulation = 6
Eigval <- c()
repeat_time = 0
repeat{
  repeat_time = repeat_time + 1
  aij <- matrix(rgamma(nn,shape = 1,scale = U[1]),
                nrow = number_of_polulation,ncol = number_of_polulation)
  bij <- matrix(rgamma(nn,shape = 1,scale = U[2]),
                nrow = number_of_polulation,ncol = number_of_polulation)
  diag(aij) <- 0
  diag(bij) <- 0
  r = runif(number_of_polulation,0.5,1)
  pars <- c(aij,bij,r)
  NN <- runif(number_of_polulation,0,1)
  tryCatch({
    out <- ode(NN, times, GLVmod, pars)
    J <- jacobian(function(x){as.vector(GLVmod(Time = 1, n = x, Pars = pars)[[1]])}, 
                  out[nrow(out),2:ncol(out)], method = "simple")
    Eigval <<- c(Eigval, eigen(J)$values[which(Re(eigen(J)$values) == max(Re(eigen(J)$values)))])
  }, error = function(e){repeat_time <<- repeat_time - 1})
  if(repeat_time >= 5000){
    break
  }
}

mydf3 <- data.frame(x = Re(Eigval), y = Im(Eigval))
myp3 <- ggplot(mydf3, aes(x = x, y = y)) +
  geom_point(color = "red", size = 1) +
  labs(
    title = "n = 6",
    x = "Real",
    y = "Imaginary"
  ) +
  theme_minimal()

myp4 <- grid.arrange(myp1 + xlim(-200, 50) + ylim(-800,800)+geom_vline(aes(xintercept = 0)), 
                     myp2 + xlim(-200, 50) + ylim(-800,800)+geom_vline(aes(xintercept = 0))+ylab(""), 
                     myp3 + xlim(-200, 50) + ylim(-800,800)+geom_vline(aes(xintercept = 0))+ylab(""), 
                     nrow = 1)
ggsave("F:/AcademyData/Computer/rcode/ecology_model/GNM3_laplasian/local_dynamics/eigenvalue_distribution/GNM.png", 
       myp4, width = 8*0.8, height = 2.5*0.8, dpi = 500)
ggsave("F:/AcademyData/Computer/rcode/ecology_model/GNM3_laplasian/local_dynamics/eigenvalue_distribution/GNM.pdf", 
       myp4, width = 8, height = 3)
mydf1$n <- 2
mydf2$n <- 3
mydf3$n <- 6
mydf <- rbind(mydf1,mydf2,mydf3)
write.csv(mydf, "F:/AcademyData/Computer/rcode/ecology_model/GNM3_laplasian/local_dynamics/eigenvalue_distribution/GLV_data.csv")
# number_of_polulation = 2
# Eigval <- c()
# repeat_time = 0
# repeat{
#   repeat_time = repeat_time + 1
#   aij <- matrix(rgamma(nn,shape = 1,scale = U[1]),
#                 nrow = number_of_polulation,ncol = number_of_polulation)
#   bij <- matrix(rgamma(nn,shape = 1,scale = U[2]),
#                 nrow = number_of_polulation,ncol = number_of_polulation)
#   diag(aij) <- 0
#   diag(bij) <- 0
#   r = runif(number_of_polulation,0.5,1)
#   pars <- c(aij,bij,r)
#   NN <- runif(number_of_polulation,0,1)
#   out <- ode(NN, times, GNMmod, pars)
#   J <- jacobian(function(x){as.vector(GNMmod(Time = 1, n = x, Pars = pars)[[1]])},
#                 out[nrow(out),2:ncol(out)], method = "simple")
#   Eigval <- c(Eigval, eigen(J)$values[which(Re(eigen(J)$values) == max(Re(eigen(J)$values)))])
#   if(any(Re(Eigval)>=0) | repeat_time >= 5000){
#     matplot(times, out[,2:ncol(out)], type = "l")
#     break
#   }
# }
model_type = "RCS"
number_of_polulation = 3
U <- c(0,4)
iterations = 1000
nn <- number_of_polulation^2
times <- seq(0, iterations, by = 0.1)
Eigval <- c()
repeat_time = 0
repeat{
  repeat_time = repeat_time + 1
  for (i in 1:3) {
    assign(paste0("k", i), -runif(1, U[1], U[2]))
  }
  aij <- matrix(c(-1, k1, 0, 
                  0, -1, k2,
                  k3, 0, -1), nrow = 3)
  r = runif(number_of_polulation,0.5,1)
  pars <- c(aij,r)
  NN <- runif(number_of_polulation,0,1)
  tryCatch({
    out <- ode(NN, times, GLVmod, pars)
    J <- jacobian(function(x){as.vector(GLVmod(Time = 1, n = x, Pars = pars)[[1]])}, 
                  out[nrow(out),2:ncol(out)], method = "simple")
    Eigval <<- c(Eigval, eigen(J)$values[which(Re(eigen(J)$values) == max(Re(eigen(J)$values)))])
  }, error = function(e){repeat_time <<- repeat_time - 1})
  if(repeat_time >= 5000){
    break
  }
}
mydf5 <- data.frame(x = Re(Eigval), y = Im(Eigval))
myp5 <- ggplot(mydf5, aes(x = x, y = y)) +
  geom_point(color = "red", size = 1) +
  labs(
    title = "RSP n = 3",
    x = "Real",
    y = "Imaginary"
  ) +
  theme_minimal()
myp5
ggsave("./eigenvalue_distribution/RSP.jpg", myp5, width = 3, height = 3, dpi = 1000)
write.csv(mydf5, file = "./eigenvalue_distribution/RSP.csv")

model_type = "GLV"
number_of_polulation = 10
U <- c(0,4)
iterations = 1000
nn <- number_of_polulation^2
times <- seq(0, iterations, by = 0.1)
Eigval <- c()
repeat_time = 0
repeat{
  repeat_time = repeat_time + 1
  aij <- - matrix(runif(9, U[1], U[2]), nrow = 3)
  diag(aij) <- -1
  r = runif(number_of_polulation,0.5,1)
  pars <- c(aij,r)
  NN <- runif(number_of_polulation,0,1)
  tryCatch({
    out <- ode(NN, times, GLVmod, pars)
    J <- jacobian(function(x){as.vector(GLVmod(Time = 1, n = x, Pars = pars)[[1]])}, 
                  out[nrow(out),2:ncol(out)], method = "simple")
    Eigval <<- c(Eigval, eigen(J)$values[which(Re(eigen(J)$values) == max(Re(eigen(J)$values)))])
  }, error = function(e){repeat_time <<- repeat_time - 1})
  if(repeat_time >= 5000){
    break
  }
}
mydf6 <- data.frame(x = Re(Eigval), y = Im(Eigval))
myp6 <- ggplot(mydf5, aes(x = x, y = y)) +
  geom_point(color = "red", size = 1) +
  labs(
    title = "GLV n = 3",
    x = "Real",
    y = "Imaginary"
  ) +
  theme_minimal()
myp6
ggsave(paste0("./eigenvalue_distribution/GLV_Uniform_", U[1], "_", U[2], ".jpg"), 
       myp5, width = 3, height = 3, dpi = 1000)
ggsave(paste0("./eigenvalue_distribution/GLV_Uniform_", U[1], "_", U[2], ".pdf"), 
       myp5, width = 3, height = 3)
mydf6$n <- 3
write.csv(mydf6, file = "./eigenvalue_distribution/GLV.csv")
