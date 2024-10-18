### This file is used to test the proportiion of different types of dynamics (stable, periodic, and chaotic) under a given condition of the model.
### Different types of dynamics are characterized by the criterion of eigenvalues and finite time Lyapunov exponent (FTLE)
### The FTLE is the criterion for chaos, the specific detection of chaos starts at line 90.

rm(list=ls())
path = ""
setwd(path)
library(deSolve)
library(ggplot2)
library(patchwork)
library(tidyr)
library(paletteer)
library(xlsx)
library(openxlsx)
library(numDeriv)
library(gg3D)
library(reticulate)
library(factoextra)
library(MASS)
Existence_rate_func <- function(x){
  Ex <- apply(x, 1, function(n){
    if(any(n<0.001)){
      return(0)
    }else{
      return(1)
    }
  })
  return(sum(Ex)/length(Ex))
}

GNMmod <- function(Time, n, Pars) {
  with( as.list(Pars), {
    n[n < 0] <- 0
    dn <- r*n*(1 - n + bij%*%n - aij%*%(n^2))
    return(list(dn))
  })
}

my_kmeans <- function(x, i){
  kmeans(x, i, nstart = 1)
}

## parameters ##
U <- c(0.25,0.25)
number_of_population = 8
iterations <- 200
nn <- number_of_population^2
times <- seq(0, iterations, by = 0.01)
last_iterations <- floor(0.75*length(times)):length(times)
precise_par <- 20

dir.create(paste("./GNM_beta", 
                 round(U[1],2), round(U[2],2), 
                 paste0("n=",number_of_population), sep = "_"))
test_time <- 0
single_equilibrium <- 0
multi_equilibrium <- 0
Osci <- 0
Chaos <- 0

mydatadf <- data.frame()

repeat{
  test_time <- test_time + 1
  set.seed(NULL)

  aij <- matrix(rgamma(nn,shape = 1,scale = U[1]),
                nrow = number_of_population,ncol = number_of_population)
  bij <- matrix(rgamma(nn,shape = 1,scale = U[2]),
                nrow = number_of_population,ncol = number_of_population)
  diag(aij) <- 0
  diag(bij) <- 0
  r = runif(number_of_population,0.5,1)
  pars <- c(aij,bij,r)
  for (i in 0:20) {
    NN <- runif(number_of_population,0,1)
    assign(paste0("out", i), 
           tryCatch({ode(NN, times, GNMmod, pars)}, 
                    error = function(e){ode(NN, times, GNMmod, pars, method = "ode45")}))
  }
  myaxessp <- order(apply(out0[last_iterations,2:ncol(out0)],2,mean), decreasing = TRUE)[1:3]+1
  out_data <- gather(data.frame(out0[,2:ncol(out0)]),key = "species", value = "density")
  out_data$time <- times

  ## stability Detection ##
  J <- jacobian(function(x){as.vector(GNMmod(Time = 1, n = x, Pars = pars)[[1]])}, 
                apply(out0[last_iterations, 2:ncol(out0)],2,mean), method = "simple")
  mymaxeigval <- max(Re(eigen(J)$values))

  ## Chaos Detection ##
  myftle_vec <- c()
  for(j in 1:20){
    out <- get(paste0("out",j))
    startings <- matrix(runif(number_of_population*precise_par, -1e-6, 1e-6), nrow = number_of_population) + 
      matrix(rep(out[ceiling(nrow(out)/2),2:ncol(out)], precise_par), nrow = number_of_population)
    endings <- matrix(0, nrow = number_of_population, ncol = precise_par)
    for (i in 1:precise_par) {
      outt <- tryCatch({ode(y = startings[,i], times = times[1:(length(times)/2)], 
                            func = GNMmod, parms = pars)}, 
                       error = function(e){ode(y = startings[,i], times = times[1:(length(times)/2)], 
                                               func = GNMmod, parms = pars, method = "ode45")}) 
      endings[,i] <- outt[nrow(outt), 2:ncol(outt)]
    }
    xbar <- apply(startings, 1, mean)
    ybar <- apply(endings, 1, mean)
    delta_x <- apply(startings, 2, function(x){x - xbar})
    delta_y <- apply(endings, 2, function(x){x - ybar})
    myJ <- delta_y %*% (ginv(t(delta_x)%*%delta_x) %*% t(delta_x))
    myftle_vec <- c(myftle_vec, log(max(svd(myJ)$d))/iterations)
  }
  max_myftle <- max(myftle_vec)
  min_myftle <- min(myftle_vec)
  
  if(mymaxeigval < 0){
    mydf <- t(data.frame(out0[nrow(out0),2:ncol(out0)]))
    for (j in 1:20) {
      out_rep <- get(paste0("out",j))
      mydf <- round(rbind(mydf, out_rep[nrow(out_rep), 2:ncol(out_rep)]),4)
    }
    best_cluster_number <- tryCatch({
      best_cluster <- fviz_nbclust(mydf, my_kmeans, k.max = 2)
      which(best_cluster$data[2] == max(best_cluster$data[2]))
    }, error = function(e){1})
  } else {
    best_cluster_number <- 0
  }
  local_dynamic <- ggplot(out_data, aes(time, density, group = species, colour = species))+
    geom_line(data = out_data[1:(iterations * 100 * number_of_population + number_of_population),], linewidth = 0.5, show.legend = FALSE) +
    geom_line(data = out_data[(iterations * 100 * number_of_population + number_of_population + 1):nrow(out_data),], show.legend = FALSE, colour = "black", linewidth = 0.5)+
    labs(x = "Time", y = "Density") +
    scale_color_manual(values = paletteer_c("grDevices::Blue-Yellow", number_of_population)) +
    theme_classic()
  path_data <- as.data.frame(out0[,myaxessp])
  colnames(path_data) <- c("sp1", "sp2", "sp3")
  mypoincare <- ggplot(path_data,aes(x=sp1, y=sp2, z=sp3)) +
    axes_3D() +
    stat_3D(geom = "path", color = "#1699AFFF") +
    labs_3D(labs = c("Sp1","Sp2","Sp3")) +
    theme_void()
  Survival_rate <- Existence_rate_func(t(out[last_iterations, 2:ncol(out)]))
  ggsave(paste(paste("./GNM_beta", 
                     round(U[1],2), round(U[2],2), 
                     paste0("n=",number_of_population), sep = "_"), 
               "/", paste(test_time, "GNM_beta", 
                          round(U[1],2), round(U[2],2), 
                          paste0("n=",number_of_population), ".png",sep = "_"), 
               sep = ""), 
         gridExtra::grid.arrange(local_dynamic,mypoincare,nrow = 1), width = 5.7, height = 3)
  mydatadf <- rbind(mydatadf, c(Survival_rate, mymaxeigval, best_cluster_number, max_myftle, min_myftle))#, myFTLE
  
  if(test_time>=500){
    break
  }
}
colnames(mydatadf) <- c("Survival Rate", "Max Eigenvalue", "Equilibrium Number", "max_FTLE", "min_FTLE")
mydatadf$dy_type <- NA
for(i in 1:nrow(mydatadf)){
  if(mydatadf$`Equilibrium Number`[i] == 1){
    mydatadf$dy_type[i] <- 1
  } else if(mydatadf$`Equilibrium Number`[i] == 2){
    mydatadf$dy_type[i] <- 2
  } else{
    mydatadf$dy_type[i] <- 0
  }
}
write.csv(mydatadf, 
          paste(paste("./GNM_beta", 
                      round(U[1],2), round(U[2],2), 
                      paste0("n=",number_of_population), sep = "_"), 
                "/", paste("GNM_beta", 
                           round(U[1],2), round(U[2],2), 
                           paste0("n=",number_of_population), "data.csv",sep = "_"), 
                sep = ""))


#   myftle <- c(myftle, log(sqrt(max(eigen(t(Delta_matrix)%*%Delta_matrix)$values)))/iterations)
#   if(repeat_time >= 20){
#     break
#   }
# }
# myFTLE <- mean(myftle)
