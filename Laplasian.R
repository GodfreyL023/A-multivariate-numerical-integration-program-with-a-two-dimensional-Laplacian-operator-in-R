rm(list=ls())
path = "Set work path here"
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

Your_Model_2D <- function(time, state, parms, N, Da, dx, dy) { 
  n <- array(state, dim = c(length(state)/(N^2), N, N))
  ### For a system with n variables dispersing in a N*N space, the initial input should be a vector transfered from n*N*N array.
  ### Then the vector as the initial value is transfered back to the array.
  with( as.list(parms), { 

    ### Design your own model instead of the next 3 rows
    positive_interaction <- apply(n, c(2,3), function(x){bij %*% x})
    netative_interaction <- apply(n, c(2,3), function(x){-aij %*% (x^2)})
    dn <- r*n*(1 - n + positive_interaction + netative_interaction)
    ### Here I exhibit a system in accordance with our article, which is a variation of the well-known Lotka-Volterra model.
    ### Because the partial derivative of the afectted variable against exerting variables change from positive to negative as exerting variables increase, we call 
    ### this model as Generalized Non-monotonous(GNM) Interaction Model.

    ### The laplacian operator is performed as follow
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

mr_gradient <- c(0, 10^(-7:-1)) ### Here set the diffusion rate gradient you need
number_of_polulation = 6 ### Here set the dimension of the variable or number of the variables you need
nn = number_of_polulation^2

### Here set the number of iterations and integration step_length
iterations = 100 
step_length = 0.01
times <- seq(0, iterations, step_length)

### Here set the spatial resolution you need
N <- 50
dx <- 1/N
dy <- 1/N

### Desigh all the other parameters you need to input to the model
U <- c(2,2)
aij <- round(matrix(rgamma(nn, 1, scale = U[1]), nrow = number_of_polulation, 
                    ncol = number_of_polulation),4)#
diag(aij)=0
bij <- round(matrix(rgamma(nn, 1, scale = U[2]), nrow = number_of_polulation, 
                    ncol = number_of_polulation),4)#
diag(bij)=0
r <- round(runif(number_of_polulation,0.5,1),4)
pars <- c(aij = aij, bij = bij, r = r)

### Here design your initial values that will be input. 
### Notice that the initial values should be stored in a number_of_polulation*N*N vector
N_ini <- runif(number_of_polulation*N*N, 0, 1)
n <- N_ini

### Traverse the diffusion rate gradient in the loop below.
### The parameters set above can also be set inside the loop, then each loop won't share the same group of parameters. 

taken_time <- 0 ### timer
for (i in 1:length(mr_gradient)) {
  start_time <- Sys.time()
  Da <- mr_gradient[i]
  out <- ode.2D(y = N_ini, times = times, func = GNMmod2D, parms = pars, 
                dimens = c(N, N), N = N, dx = dx, dy = dy, Da = Da, ynames = FALSE, hini = 0.01, method = "ode45")
  ### Change and adjust the integrator according the deSolve profile, such as change the method to "ode23","lsoda". 

  ### Transfer the final output of the integrator to the array in accordance with the initial designation
  premydata <- array(as.vector(out[,2:ncol(out)]), dim = c(iterations+1, number_of_polulation,N^2))
  mydata <- aperm(premydata,c(2,3,1))

  ### Figures exhibit the dynamics and the final spatial pattern
  mydata1 <- round(apply(mydata,2:3,sum),4)
  data <- gather(data.frame(t(mydata1[,1:(iterations+1)])),key = "grid",value = "abundance")
  data$times <- rep(1:(iterations+1),N^2)
  p1 <- ggplot(data,aes(times,abundance),group = grid)+
    geom_line(aes(color = grid),show.legend = FALSE,linewidth = 0.7) +
    labs(x = "Time", y = "abundance") +
    theme_classic()
  
  data <- data.frame(mydata1[,iterations+1],1:N,rep(1:N,each=N))
  data <- round(data,4)
  names(data) <- c("abundance","x","y")
  data$abundance <- (data$abundance)/(max(data$abundance))
  p2 <- ggplot(data,aes(x,y,fill = abundance))+
    geom_tile(color = "white", size = 0.25) +
    scale_fill_gradient2(low = "#00BFC4",mid = "white",high = "#F8766D",midpoint = 0.5,limits = c(0,1),name = "",,na.value = "#F8766D") +
    labs(x = "", y = "") +
    theme_classic()
  p3 <- grid.arrange(p1,p2,nrow=2)
  ggsave(paste(number_of_polulation, "Gamma", U[1], U[2],
               "mr=",Da,"Community.png",sep ="_"),
         p3,dpi=1000,width = 4,height = 6.4)
  
  #population level figs
  interest_sp <- apply(mydata[,,(iterations*0.9):(iterations+1)],1,mean)
  interest_sp <- which(interest_sp == max(interest_sp))
  interest_sp <- interest_sp[1]
  mydata2 <- round(mydata[interest_sp,,],4)
  data <- gather(data.frame(t(mydata2[,1:(iterations+1)])),key = "grid",value = "abundance")
  data$times <- rep(1:(iterations+1),N^2)
  p4 <- ggplot(data,aes(times,abundance),group = grid)+
    geom_line(aes(color = grid),show.legend = FALSE,linewidth = 0.7) +
    labs(x = "Time", y = "abundance") +
    theme_classic()
  
  data <- data.frame(mydata2[,iterations+1],1:N,rep(1:N,each=N))
  data <- round(data,4)
  names(data) <- c("abundance","x","y")
  data$abundance <- (data$abundance)/(max(data$abundance))
  p5 <- ggplot(data,aes(x,y,fill = abundance))+
    geom_tile(color = "white", size = 0.25) +
    scale_fill_gradient2(low = "#00BFC4",mid = "white",high = "#F8766D",midpoint = 0.5,limits = c(0,1),name = "",na.value = "#F8766D") +
    labs(x = "", y = "") +
    theme_classic()
  
  p6 <- grid.arrange(p4,p5,nrow=2)
  ggsave(paste(number_of_polulation, "Gamma", U[1], U[2],
               "mr=",Da,"Population.png",sep ="_"),
         p6,dpi=1000,width = 4,height = 6.4)
  
  #Gini index that shows the evenness of the final spatial pattern
  gini <- function(x) {
    G <- 2 * sum((1:length(x)) * sort(x)) / length(x) - 1 - 1 / length(x)
    G[G<0]=0
    return(G)
  }
  
  xC <- apply(mydata[,,],c(2,3),sum)
  xxC <- apply(xC,2,sort)
  xxxC <- tryCatch({
    apply(xxC, 2, prop.table)
  },
  error = function(e){
    xxC <- do.call(cbind,xxC)
    apply(xxC, 2, prop.table)
  })
  G_Community <- apply(xxxC, 2, gini)
  
  xP <- apply(mydata[interest_sp,,],2,sort)
  xxP <- tryCatch({
    apply(xP, 2, prop.table)
  },
  error = function(e){
    xP <- do.call(cbind,xP)
    apply(xP, 2, prop.table)
  })
  
  G_Population <- apply(xxP, 2, gini)
  length1 <- length(G_Community)
  length2 <- length(G_Population)
  max_length <- max(length1, length2)
  if (length1 < max_length) {
    G_Community <- c(G_Community, rep(NA, max_length - length1))
  } else if (length2 < max_length) {
    G_Population <- c(G_Population, rep(NA, max_length - length2))
  }
  
  mydata3 <- tryCatch({
    round(data.frame(G_Community,G_Population),4)
  },
  error = function(e){
    GG <- t(G_Population)
    round(data.frame(G_Community,GG[interest_sp,]),4)
  })
  
  names(mydata3) <- c("Community","Population")
  data1 <- gather(mydata3[1:(iterations+1),],key = "Level",value="Value")
  data1$Time <- 1:(iterations+1)
  
  #Moran's I shows the spatial autocorrelation of the space
  o <- seq(1,N^2,by=1)
  w_matrix <- matrix(0,nrow = N^2, ncol = N^2)
  
  for (i in 1:N^2) { # tranverse grids
    x = ceiling(o[i]/N) #x-axis in matrix
    y = o[i] %% N #y-axis in matrix
    if(y == 0){
      y = N
    }
    
    ddd = 1:N
    
    if((y-1) %in% ddd == TRUE){
      w_matrix[i-1,i] = 1
    }
    
    if((y+1) %in% ddd == TRUE){
      w_matrix[i+1,i] = 1
    }
    
    if((x-1) %in% ddd == TRUE){
      w_matrix[i,i-N] = 1
    }
    
    if((x+1) %in% ddd == TRUE){
      w_matrix[i,i+N] = 1
    }
  }
  
  getMoranI <- function(x){
    tryCatch({
      MI <- Moran.I(x,w_matrix)
      return(unlist(MI))
    },
    error = function(e){
      EE <- -1/(length(x)-1)
      return(c(1,EE,0,0))
    })
  }
  
  my_I_P <- round(apply(mydata[interest_sp,,],2,getMoranI),4)
  my_I_P_mean <- apply(my_I_P[,(iterations*0.9):(iterations+1)], 1, mean)
  
  my_I_C <- round(apply(apply(mydata[,,],2:3, sum), 2, getMoranI),4)
  my_I_C_mean <- apply(my_I_C[,(iterations*0.9):(iterations+1)], 1, mean)
  
  IC_paste <- round(my_I_C_mean[4],3)
  if(IC_paste<0.01){
    IC_paste <- paste("< 0.01 **")
  }else if(IC_paste<0.05){
    IC_paste <- paste("=",IC_paste,"*")
  }else{
    IC_paste <- paste("=",IC_paste)
  }
  IP_paste <- round(my_I_P_mean[4],3)
  if(IP_paste<0.01){
    IP_paste <- paste("< 0.01 **")
  }else if(IP_paste<0.05){
    IP_paste <- paste("=",IP_paste,"*")
  }else{
    IP_paste <- paste(IP_paste)
  }
  
  mydata4 <- data.frame(my_I_P[1,],my_I_C[1,])
  names(mydata4)<-c("Population","Community")
  data2 <- gather(mydata4[1:(iterations+1),],key="Level",value = "Value")
  data2$Time <- 1:(iterations+1)
  
  
  data <- rbind(data1,data2)
  data$Index <- c(rep("Gini Index", nrow(data1)), rep("Moran's I",nrow(data2)))
  
  text_df <- data.frame(x = c(0,0),
                        y=c(max(data$Value)*1.05, max(data$Value)*1.1),
                        Index=c("Gini Index", "Moran's I"),
                        label=c("", paste("P Value of Moran's I", IC_paste,",", IP_paste)))
  p10 <- ggplot(data,aes(Time,Value,group=Level,color=Level))+
    geom_line(linewidth = 0.7)+
    facet_wrap(~Index,dir = "v",strip.position="right",scales = "free_y")+
    geom_text(data=text_df, aes(x=x,y=y,label=label),
              color="black", hjust=0, inherit.aes = FALSE)+
    theme_classic()
  
  ggsave(paste(N,number_of_polulation, number_of_polulation, "Gamma", U[1], U[2],
               "mr=", Da,"Gini and Moran's I.png", sep ="_"),
         p10,dpi=1000,width = 5,height = 5)
  
  pvaluedf <- data.frame(c(round(my_I_C_mean[4],3)),round(my_I_P_mean[4],3))
  names(pvaluedf) <- c("Community","Population")
  
  #Data archiving
  mybook <- createWorkbook()
  C_data <- addWorksheet(mybook,sheetName = "Community Data")
  P_data <- addWorksheet(mybook,sheetName = "Population Data")
  aij_r <- addWorksheet(mybook,sheetName = "aij and r")
  ini_n <- addWorksheet(mybook,sheetName = "Initial Populaiton Sizes")
  gini <- addWorksheet(mybook,sheetName = "Gini Index")
  MoranI <- addWorksheet(mybook,sheetName = "Moran's I")
  P_moran <- addWorksheet(mybook,sheetName = "P of Moran's I")
  
  writeData(mybook,sheet = C_data,mydata1)
  writeData(mybook,sheet = P_data,mydata2)
  writeData(mybook,sheet = aij_r,data.frame(r,aij))
  writeData(mybook,sheet = ini_n,data.frame(t(mydata[,,1])))
  writeData(mybook,sheet = gini,mydata3)
  writeData(mybook,sheet = MoranI,mydata4)
  writeData(mybook,sheet = P_moran,pvaluedf)
  
  saveWorkbook(mybook, paste(N, number_of_polulation, "Gamma", U[1], U[2],
                             "mr=", Da, ".csv", sep ="_"))
  end_time <- Sys.time()
  taken_time <- taken_time + end_time - start_time
  print(taken_time)
 }






