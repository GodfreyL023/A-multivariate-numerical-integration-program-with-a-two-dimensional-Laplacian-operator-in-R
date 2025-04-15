rm(list=ls())
path = "F:/AcademyData/Computer/rcode/ecology_model/Stochasticity/2sp_GNM_test4/"
tryCatch({
  dir.create(path)
}, error = function(e){})
setwd(path)
library(deSolve)
library(ggplot2)
library(gridExtra)
library(abind)
library(tidyr)
library(xlsx)
library(openxlsx)
library(ape)
library(paletteer)

GNMmod2D <- function(time, state, parms, N, Da, dx, dy) { 
  state[state<=0] <- 0
  n <- array(state, dim = c(length(state)/(N^2), N, N))
  with( as.list(parms), { 
    
    dn <- n * (1 + rnorm(length(state)/(N^2), 0, sigma)) * (1 - n + 
                 apply(n, c(2,3), function(x){bij %*% x}) - 
                 apply(n, c(2,3), function(x){aij %*% x^2}))
    
    zero <- numeric(N)
    Flux_x <- abind(matrix(rep(zero,dim(n)[1]), nrow = dim(n)[1], byrow = TRUE), 
                    (-Da * (n[,-1,] - n[,-N,]) / dx), 
                    matrix(rep(zero,dim(n)[1]), nrow = dim(n)[1], byrow = TRUE), 
                    along = 2)
    Flux_y <- abind(matrix(rep(zero,dim(n)[1]), nrow = dim(n)[1], byrow = TRUE), 
                    -Da * (n[,,-1] - n[,,-N]) / dy, 
                    matrix(rep(zero,dim(n)[1]), nrow = dim(n)[1], byrow = TRUE), 
                    along = 3)
    dn <- dn - (Flux_x[,-1,] - Flux_x[,-(N+1),]) / dx - 
      (Flux_y[,,-1] - Flux_y[,,-(N+1)]) / dy
    
    return(list(as.vector(dn)))
  })
}

mr_gradient <- c(0, sort(as.vector(outer(c(1,3,5,7), 10^(-6:-3)))), 0.01, 0.1)# 
number_of_population = 2
nn = number_of_population^2
iterations = 200
step_length <- 0.1
orbits_proportion <- 0.1

a1 <- 2
b1 <- 2
aij <- matrix(c(0, a1, b1, 0), nrow = 2, byrow = TRUE)
bij <- aij
sigma <- 0
N <- 5
dx <- 1/N
dy <- 1/N


times <- seq(0, iterations, step_length)
pars <- c(aij, bij, sigma)
N_ini <- runif(number_of_population*N*N, 0.01, 1)
taken_time <- 0
for (i in 1:length(mr_gradient)) {
  start_time <- Sys.time()
  Da <- mr_gradient[i]
  out <- ode.2D(y = N_ini, times = times, func = GNMmod2D, parms = pars, 
                dimens = c(N, N), N = N, dx = dx, dy = dy, Da = Da, 
                ynames = FALSE, method = "ode45", 
                atol = 0.001, rtol = 0.001)
  
  premydata <- array(as.vector(out[,2:ncol(out)]), dim = c(iterations/step_length + 1, number_of_population, N^2))
  mydata <- aperm(premydata,c(2,3,1))
  rm(premydata)
  rm(out)
  gc()
  # matplot(premydata[,6,],type = "l")
  # matplot(t(mydata[6,,]),type = "l")
  
  orbits_number <- sample(1:(N^2), ceiling(orbits_proportion * N^2))
  if(length(orbits_number)<2){
    orbits_number <- sample(1:(N^2), ceiling(orbits_proportion * N^2)+1)
  }
  mydata1 <- apply(mydata,2:3,sum)
  data <- gather(data.frame(t(mydata1[orbits_number,1:(iterations/step_length+1)])),key = "grid",value = "abundance")
  data$times <- times
  p1 <- ggplot(data, aes(times,abundance),group = grid)+
    geom_line(aes(color = grid),show.legend = FALSE,linewidth = 0.7) +
    scale_x_log10(limits = c(1,iterations))+
    #scale_y_log10()+
    labs(x = "Time", y = "abundance") +
    theme_classic()
  
  data <- data.frame(mydata1[,iterations/step_length+1],1:N,rep(1:N,each=N))
  data <- round(data,4)
  names(data) <- c("abundance","x","y")
  data$abundance <- (data$abundance)/(max(data$abundance))
  p2 <- ggplot(data,aes(x,y,fill = abundance))+
    geom_tile(color = "white", size = 0.25) +
    scale_fill_gradient2(low = "#2D3184FF",mid = "#1296AEFF",high = "yellow",midpoint = 0.5,limits = c(0,1),name = "",,na.value = "yellow") +
    labs(x = "", y = "") +
    theme_classic()
  p3 <- grid.arrange(p1,p2,nrow=2)
  ggsave(paste(number_of_population, "Noise", a1, b1,
               "mr=",Da,"sigma=",sigma,"Community.png",sep ="_"),
         p3,dpi=300,width = 4,height = 6.4)
  
  #population level figs
  interest_sp <- apply(mydata[,,(iterations*0.9/step_length):(iterations/step_length+1)],1,mean)
  interest_sp1 <- order(unique(interest_sp), decreasing = TRUE)[1]
  interest_sp2 <- order(unique(interest_sp), decreasing = TRUE)[2]
  mydata2 <- mydata[interest_sp1,,]
  data <- gather(data.frame(t(mydata2[orbits_number,1:(iterations/step_length+1)])),key = "grid",value = "abundance")
  data$times <- times
  p4 <- ggplot(data,aes(times,abundance),group = grid)+
    geom_line(aes(color = grid),show.legend = FALSE,linewidth = 0.7) +
    scale_x_log10(limits = c(1,iterations))+
    #scale_y_log10()+
    labs(x = "Time", y = "abundance") +
    theme_classic()
  
  data <- data.frame(mydata2[,iterations/step_length+1],1:N,rep(1:N,each=N))
  data <- round(data,4)
  names(data) <- c("abundance","x","y")
  data$abundance <- (data$abundance)/(max(data$abundance))
  p5 <- ggplot(data,aes(x,y,fill = abundance))+
    geom_tile(color = "white", size = 0.25) +
    scale_fill_gradient2(low = "#2D3184FF",mid = "#1296AEFF",high = "yellow",midpoint = 0.5,limits = c(0,1),name = "",na.value = "yellow") +
    labs(x = "", y = "") +
    theme_classic()
  
  p6 <- grid.arrange(p4,p5,nrow=2)
  ggsave(paste(number_of_population, "Noise", a1, b1, 
               "mr=",Da, "sigma=",sigma, "Population1.png",sep ="_"),
         p6,dpi=300,width = 4,height = 6.4)
  
  #Data archiving
  mybook <- createWorkbook()
  C_data <- addWorksheet(mybook,sheetName = "Community Data")
  P_data <- addWorksheet(mybook,sheetName = "Population Data")
  aij_bij <- addWorksheet(mybook,sheetName = "aij and bij")
  ini_n <- addWorksheet(mybook,sheetName = "Initiald Populaiton Sizes")
  chunk_size <- 100
  mydata11 <- t(round(mydata1[,seq(1,ncol(mydata1),by = 5)],4))
  rm(mydata1)
  gc()
  for (j in seq(1, nrow(mydata11), by = chunk_size)) {
    data <- mydata11[j:min(j+chunk_size-1, nrow(mydata11)),]
    if(is.vector(data)){
      data <- matrix(data, nrow = 1)
    }
    writeData(mybook, sheet = C_data, data, startRow = j, colNames = FALSE)
  }
  mydata22 <- t(round(mydata2[,seq(1,ncol(mydata2),by = 5)],4))
  rm(mydata2)
  gc()
  for (j in seq(1, nrow(mydata22), by = chunk_size)) {
    data <- mydata22[j:min(j+chunk_size-1, nrow(mydata11)),]
    if(is.vector(data)){
      data <- matrix(data, nrow = 1)
    }
    writeData(mybook, sheet = P_data, data, startRow = j, colNames = FALSE)
  }
  #writeData(mybook,sheet = P_data,round(mydata2[,seq(1,ncol(mydata2),by = 10)],4))
  writeData(mybook,sheet = aij_bij,rbind(data.frame(aij), data.frame(bij)))
  writeData(mybook,sheet = ini_n,data.frame(t(mydata[,,1])))
  
  saveWorkbook(mybook, paste(N, number_of_population, "2sp_GNMsto", a1, b1,
                             "mr=", Da, "sigma=", sigma, ".xlsx", sep ="_"))
  end_time <- Sys.time()
  taken_time <- taken_time + end_time - start_time
  print(taken_time)
}
