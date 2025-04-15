rm(list=ls())
path = "F:/AcademyData/Computer/rcode/ecology_model/GNM3_laplasian/preliminary_experiment/GNMRSP_test33/"
tryCatch({
  dir.create(path)
}, error = function(e){})
setwd(path)
library(deSolve)
library(ggplot2)
library(gridExtra)
library(abind)
library(tidyr)
library(dplyr)
library(xlsx)
library(openxlsx)
library(ape)
library(paletteer)

GNMmod2D <- function(time, state, parms, N, Da, dx, dy) { 
  state[which(state<0)] <- 0
  n <- array(state, dim = c(length(state)/(N^2), N, N))
  with( as.list(parms), { 
    
    dn <- n * (1 - n + 
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

mr_gradient <- c(3e-5, 5e-5, 7e-5, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4)# 0, sort(as.vector(outer(c(1,5), 10^(-6:-3)))), 0.01, 0.1
number_of_population = 3
nn = number_of_population^2
iterations = 500
step_length <- 0.1
orbits_proportion <- 0.05

for (i in 1:3) {
  assign(paste0("a", i), c(5,0.1,3)[i])
}
aij <- matrix(c(0, a1, a2,
                a2, 0, a1,
                a1, a2, 0), nrow = 3, byrow = TRUE)
bij <- aij
N <- 50
dx <- 1/N
dy <- 1/N

times <- seq(0, iterations, step_length)
pars <- c(aij, bij)
N_ini <- runif(number_of_population*N*N, 0, 1)
taken_time <- 0
for (i in 1:length(mr_gradient)) {
  start_time <- Sys.time()
  Da <- mr_gradient[i]
  out <- ode.2D(y = N_ini, times = times, func = GNMmod2D, parms = pars, 
                dimens = c(N, N), N = N, dx = dx, dy = dy, Da = Da, ynames = FALSE, method = "ode45", 
                atol = 0.0001, rtol = 0.0001)#, hini = 0.01
  
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
  orbit_number <- sort(orbits_number)
  mydata1 <- apply(mydata,2:3,sum)
  data <- gather(data.frame(t(mydata1[orbits_number,1:(iterations/step_length+1)])),key = "grid",value = "abundance")
  data$times <- times
  p1 <- ggplot(data, aes(times,abundance),group = grid)+
    geom_line(aes(color = grid),show.legend = FALSE,linewidth = 0.7) +
    #scale_x_log10()+
    #scale_y_log10()+
    labs(x = "Time", y = "Density") +
    theme_classic()
  
  data <- data.frame(mydata1[,iterations/step_length+1],1:N,rep(1:N,each=N))
  data <- round(data,4)
  names(data) <- c("abundance","x","y")
  data$abundance <- (data$abundance)/(max(data$abundance))
  data$x <- data$x / N
  data$y <- data$y / N
  p2 <- ggplot(data,aes(x,y,fill = abundance))+
    geom_tile(color = "white", size = 0.25) +
    scale_fill_gradient2(low = "#2D3184FF",mid = "#1296AEFF",high = "yellow",midpoint = 0.5,limits = c(0,1),name = "",,na.value = "yellow") +
    labs(x = "", y = "") +
    theme_classic()
  p3 <- grid.arrange(p1,p2,nrow=2)
  ggsave(paste(number_of_population, "Determine", a1, a2, a3,
               "mr=",Da,"Community.png",sep ="_"),
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
    #scale_x_log10()+
    #scale_y_log10()+
    labs(x = "Time", y = "Density") +
    theme_classic()
  
  data <- data.frame(mydata2[,iterations/step_length+1],1:N,rep(1:N,each=N))
  data <- round(data,4)
  names(data) <- c("abundance","x","y")
  data$abundance <- (data$abundance)/(max(data$abundance))
  data$x <- data$x / N
  data$y <- data$y / N
  p5 <- ggplot(data,aes(x,y,fill = abundance))+
    geom_tile(color = "white", size = 0.25) +
    scale_fill_gradient2(low = "#2D3184FF",mid = "#1296AEFF",high = "yellow",midpoint = 0.5,limits = c(0,1),name = "",na.value = "yellow") +
    labs(x = "", y = "") +
    theme_classic()
  
  p6 <- grid.arrange(p4,p5,nrow=2)
  ggsave(paste(number_of_population, "Determine", a1, a2, a3,
               "mr=",Da,"Population1.png",sep ="_"),
         p6,dpi=300,width = 4,height = 6.4)
  
  dataa <- data.frame()
  mydata3 <- mydata[interest_sp2,,]
  for (i in 1:length(orbits_number)) {
    dataa <- rbind(dataa, data.frame(x = mydata2[orbits_number[i],],
                                     y = mydata3[orbits_number[i],], 
                                     space = i))
  }
  
  arrows <- c()
  for (i in 1:length(orbits_number)) {
    arrows <- c(arrows, 2 + (i-1)*(iterations/step_length + 1))
  }
  arrow_data <- dataa[arrows,1:2]
  colnames(arrow_data) <- c("x_start", "y_start")
  arrow_data$x_end <- dataa[arrows + 1,1]
  arrow_data$y_end <- dataa[arrows + 1,2]
  arrow_data$space <- 1:nrow(arrow_data)
  
  p7 <- ggplot() +
    geom_path(data = dataa, aes(x, y, group = space, colour = space))+  
    geom_segment(
      data = arrow_data,
      aes(x = x_start, y = y_start, xend = x_end, yend = y_end, colour = space),
      arrow = arrow(
        type = "closed",          # 实心箭头
        length = unit(0.1, "cm")  # 箭头长度
      ),           # 箭头宽度
    )

  ggsave( paste(number_of_population, "Determine", a1, a2, a3, 
                "mr=",Da,"Population_Path.png",sep ="_"),p7 , width = 4, height = 3.2)
  
  #Data archiving
  mybook <- createWorkbook()
  C_data <- addWorksheet(mybook,sheetName = "Community Data")
  P1_data <- addWorksheet(mybook,sheetName = "Population1 Data")
  P2_data <- addWorksheet(mybook,sheetName = "Population2 Data")
  aij_bij <- addWorksheet(mybook,sheetName = "aij and bij")
  ini_n <- addWorksheet(mybook,sheetName = "Initiald Populaiton Sizes")
  chunk_size <- 10
  #rm(mydata1)
  gc()
  for (j in seq(1, nrow(mydata1), by = chunk_size)) {
    data <- mydata1[j:min(j+chunk_size-1, nrow(mydata1)),]
    if(is.vector(data)){
      data <- matrix(data, nrow = 1)
    }
    writeData(mybook, sheet = C_data, data, startRow = j, colNames = FALSE)
  }
  #rm(mydata2)
  gc()
  for (j in seq(1, nrow(mydata2), by = chunk_size)) {
    data <- mydata2[j:min(j+chunk_size-1, nrow(mydata2)),]
    if(is.vector(data)){
      data <- matrix(data, nrow = 1)
    }
    writeData(mybook, sheet = P1_data, data, startRow = j, colNames = FALSE)
  }
  gc()
  for (j in seq(1, nrow(mydata3), by = chunk_size)) {
    data <- mydata3[j:min(j+chunk_size-1, nrow(mydata3)),]
    if(is.vector(data)){
      data <- matrix(data, nrow = 1)
    }
    writeData(mybook, sheet = P2_data, data, startRow = j, colNames = FALSE)
  }

  #writeData(mybook,sheet = P_data,round(mydata2[,seq(1,ncol(mydata2),by = 10)],4))
  writeData(mybook,sheet = aij_bij,rbind(data.frame(aij), data.frame(bij)))
  writeData(mybook,sheet = ini_n,data.frame(t(mydata[,,1])))

  saveWorkbook(mybook, paste(N, number_of_population, "GNMRSP", a1, a2, a3,
                             "mr=", Da, ".xlsx", sep ="_"))
  end_time <- Sys.time()
  taken_time <- taken_time + end_time - start_time
  print(taken_time)
}






