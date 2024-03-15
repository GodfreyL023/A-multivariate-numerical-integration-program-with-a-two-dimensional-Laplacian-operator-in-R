### Use this program to manipulate massive simulations and drag out the general features of your model.

### Since R is born not good at intense computation, I highly recommend run this program in Intel MKL environment or after applying any other 
### multi-core arithmetic acceleration method. And even though the programme will be really time-consuming because the system is too complex, i.e., 
### too many dimensions. Therefore substantial preliminary experiments and stable computing equipments are necessary.

rm(list=ls())

setwd("Set your own work direction")
library(deSolve)
library(ggplot2)
library(gridExtra)
library(abind)
library(tidyr)
library(xlsx)
library(openxlsx)
library(progress)
options(scipen = 100) ### Make R do not use scientific notation to read and write files correctly.

### The model is the same as "Laplacian.R", please check the file before using this program.
GNMmod2D <- function(time, state, parms, N, Da, dx, dy) { 
  n <- array(state, dim = c(length(state)/(N^2), N, N))
  n[n < 0] <- 0
  # n <- array(data, dim = c(num, x, y))
  with( as.list(parms), { 
    
    positive_interaction <- apply(n, c(2,3), function(x){bij %*% x})
    netative_interaction <- apply(n, c(2,3), function(x){-aij %*% (x^2)})
    dn <- r*n*(1 - n + positive_interaction + netative_interaction)
    
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

### Here we define the function to calculate the diversity or the survival rate of each space
### You can easily define your own function to detect the features you are interested in since the output data is still in the form of n*N*N array.
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


mr_gradient <- c(10^(-4)) ###
target_terms <- c("Biomass","Existence","Osci","Biomass_var",
                  "Existence_var","Osci_var","Aij", "Bij","r","N_ini")#, "raw_data"
number_of_polulation = 8
nn = number_of_polulation^2
iterations = 100
last_iterations = (0.75*iterations+1):iterations
N <- 30
dx <- 1/N
dy <- 1/N
times <- seq(0, iterations, 1)
repeat_required <- 50

par_df <- expand.grid(seq(0.25,2,by = 0.25),seq(0.25,2,by = 0.25))

for (iii in 1:length(target_terms)) {
  assign(paste0(target_terms[iii],"_xlsx"),createWorkbook())
  for (jjj in 1:length(mr_gradient)) {
    adw <- addWorksheet(get(paste0(target_terms[iii],"_xlsx")),
                        sheetName = paste0("mr=",log10(mr_gradient[jjj])))
    assign(paste0(target_terms[iii],"_mr=",log10(mr_gradient[jjj])),adw)
  }
}

total_iterations <- length(mr_gradient)*nrow(par_df)*repeat_required
pb <- progress_bar$new(
  format = "[:bar] :percent, ETA: :eta", 
  total = total_iterations  # 设置总迭代次数
)

for (j in 1:length(mr_gradient)) {
  Da <- mr_gradient[j]
  Existence_matrix <- matrix(0, nrow = nrow(par_df), ncol = repeat_required)
  row.names(Existence_matrix) <- paste0(expression("\u03B21 = "),par_df[,1],expression(" \u03B22 = "),par_df[,2])
  Osci_matrix <- matrix(0, nrow = nrow(par_df), ncol = repeat_required)
  row.names(Osci_matrix) <- paste0(expression("\u03B21 = "),par_df[,1],expression(" \u03B22 = "),par_df[,2])
  Biomass_matrix <- matrix(0, nrow = nrow(par_df), ncol = repeat_required)
  row.names(Biomass_matrix) <- paste0(expression("\u03B21 = "),par_df[,1],expression(" \u03B22 = "),par_df[,2])
  Biomass_var_matrix <- matrix(0, nrow = nrow(par_df), ncol = repeat_required)
  row.names(Biomass_var_matrix) <- paste0(expression("\u03B21 = "),par_df[,1],expression(" \u03B22 = "),par_df[,2])
  Existence_var_matrix <- matrix(0, nrow = nrow(par_df), ncol = repeat_required)
  row.names(Existence_var_matrix) <- paste0(expression("\u03B21 = "),par_df[,1],expression(" \u03B22 = "),par_df[,2])
  Osci_var_matrix <- matrix(0, nrow = nrow(par_df), ncol = repeat_required)
  row.names(Osci_var_matrix) <- paste0(expression("\u03B21 = "),par_df[,1],expression(" \u03B22 = "),par_df[,2])
  Aij_df <- data.frame()
  Bij_df <- data.frame()
  r_df <- data.frame()
  N_ini_df <- data.frame()
  #rawdata_df <- data.frame()
  for (i in 1:nrow(par_df)) {
    U <- c(par_df[i,1], par_df[i,2])
    repeat_times <- 0
    region_existence <- c()
    Osci <- c()
    Biomass <- c()
    Biomass_var <- c()
    Existence_var <- c()
    Osci_var <- c()
    # start_time <- Sys.time()
    repeat{
      repeat_times <- repeat_times + 1
      tryCatch({
        aij <- round(matrix(rgamma(nn, 1, scale = U[1]), nrow = number_of_polulation),4)
        diag(aij)=0
        bij <- round(matrix(rgamma(nn, 1, scale = U[2]), nrow = number_of_polulation),4)
        diag(bij)=0
        r <- round(runif(number_of_polulation,0.5,1),4)
        pars <- c(aij = aij, bij = bij, r = r)
        N_ini <- runif(number_of_polulation*N*N, 0, 1)
        out <- ode.2D(y = N_ini, times = times, func = GNMmod2D, parms = pars, dimens = c(N, N),
                      N = N, dx = dx, dy = dy, Da = Da, ynames = FALSE, method = "ode45", hini = 0.1, maxsteps = 100000)
        
        premydata <- array(as.vector(out[,2:ncol(out)]), dim = c(iterations+1, number_of_polulation,N^2))
        mydata <- aperm(premydata,c(2,3,1))
        
        #productivity detector
        my_community <- apply(mydata, 2:3, sum)
        mean_total_biomass <- mean(my_community[,last_iterations])
        var_total_biomass <- sd(as.vector(my_community[,last_iterations]))
        Biomass <<- c(Biomass, mean_total_biomass)
        Biomass_var <<- c(Biomass_var, var_total_biomass)
        
        #existence detector
        existence <- apply(mydata[,,last_iterations], 2, Existence_rate_func)
        region_existence <<- c(region_existence, mean(existence))
        Existence_var <<- c(Existence_var, sd(existence))
        
        # #Diversity detector 1 -- Shannon-Wienner
        # p <- apply(my_community[,last_iterations], 2, prop.table)
        # swdi <- apply(p,2,SW)
        # swdi[is.na(swdi)] <- 0
        # swdi_mean <- sum(swdi)/length(swdi)
        # SWDI <<- c(SWDI, swdi_mean)
        # 
        # #Diversity detector 2 -- Simpson
        # simpson <- mean(apply(p, 2, function(x){1-sum(x^2)}))
        # Simpson <<- c(Simpson, simpson)
        # 
        #Oscillation detector
        osci <- apply(my_community[, last_iterations], 1, function(x){mean(abs(x-sum(x)/length(x)))})
        Osci <<- c(Osci, sum(osci)/length(osci))
        Osci_var <<- c(Osci_var,sd(osci))
        
        #Parameters archiving
        adf <- data.frame()
        adf[1,1:(nn+1)] <- c(paste0(expression("\u03B21 = "),
                             par_df[i,1],expression(" \u03B22 = "),par_df[i,2]),
                      as.vector(aij))
        Aij_df <<- rbind(Aij_df, adf)
        
        bdf <- data.frame()
        bdf[1,1:(nn+1)] <- c(paste0(expression("\u03B21 = "),
                                    par_df[i,1],expression(" \u03B22 = "),par_df[i,2]),
                             as.vector(bij))
        Bij_df <<- rbind(Bij_df, bdf)
        
        rdf <- data.frame()
        rdf[1,1:(number_of_polulation+1)] <- c(paste0(expression("\u03B21 = "),
                                    par_df[i,1],expression(" \u03B22 = "),par_df[i,2]),r)
        r_df <<- rbind(r_df, rdf)
        
        ndf <- data.frame()
        ndf[1,1:(length(N_ini)+1)] <- c(paste0(expression("\u03B21 = "),
                                                      par_df[i,1],expression(" \u03B22 = "),par_df[i,2]),round(N_ini,4))
        N_ini_df <<- rbind(N_ini_df,ndf)
        
        # rawdf <- data.frame()
        # mydata <- round(mydata,4)
        # for (sp_pool in 1:number_of_polulation) {
        #   for (grids in 1:(N^2)) {
        #     rawdf <- rbind(rawdf, c(paste0(expression("\u03B21 = "),
        #                                    par_df[i,1],expression(" \u03B22 = "),par_df[i,2]), 
        #                             repeat_times, sp_pool, grids, mydata[sp_pool,grids,]))
        #   }
        # }
        # colnames(rawdf) <- colnames(rawdata_df)
        # rawdata_df <<- rbind(rawdata_df,rawdf)
        pb$tick()
        if(repeat_times >= repeat_required){break}
      },
      error = function(e){
        repeat_times <<- repeat_times - 1
      })
    }
    Biomass_matrix[i,] <- round(Biomass,4)
    Biomass_var_matrix[i,] <- round(Biomass_var, 4)
    Existence_matrix[i,] <- round(region_existence,4)
    Existence_var_matrix[i,] <- round(Existence_var,4)
    Osci_matrix[i,] <- round(Osci,4)
    Osci_var_matrix[i,] <- round(Osci_var,4)
    
    # end_time <- Sys.time()
    # print(end_time - start_time)
  }
  writeData(Biomass_xlsx,sheet = get(paste0(target_terms[1],"_mr=",log10(Da))),
            Biomass_matrix,colNames = FALSE, rowNames = TRUE)
  writeData(Existence_xlsx,sheet = get(paste0(target_terms[2],"_mr=",log10(Da))),
            Existence_matrix,colNames = FALSE, rowNames = TRUE)
  writeData(Osci_xlsx,sheet = get(paste0(target_terms[3],"_mr=",log10(Da))),
            Osci_matrix,colNames = FALSE, rowNames = TRUE)
  writeData(Biomass_var_xlsx, sheet = get(paste0(target_terms[4],"_mr=",log10(Da))),
            Biomass_var_matrix,colNames = FALSE, rowNames = TRUE)
  writeData(Existence_var_xlsx, sheet = get(paste0(target_terms[5],"_mr=",log10(Da))),
            Existence_var_matrix,colNames = FALSE, rowNames = TRUE)
  writeData(Osci_var_xlsx, sheet = get(paste0(target_terms[6],"_mr=",log10(Da))),
            Osci_var_matrix,colNames = FALSE, rowNames = TRUE)
  writeData(Aij_xlsx, sheet = get(paste0(target_terms[7],"_mr=",log10(Da))),
            Aij_df, colNames = FALSE, rowNames = FALSE)
  writeData(Bij_xlsx, sheet = get(paste0(target_terms[8],"_mr=",log10(Da))),
            Bij_df, colNames = FALSE, rowNames = FALSE)
  writeData(r_xlsx, sheet = get(paste0(target_terms[9],"_mr=",log10(Da))),
            r_df, colNames = FALSE, rowNames = FALSE)
  writeData(N_ini_xlsx, sheet = get(paste0(target_terms[10],"_mr=",log10(Da))),
            N_ini_df, colNames = FALSE, rowNames = FALSE)
  # writeData(raw_data_xlsx, sheet = get(paste0(target_terms[11],"_mr=",log10(Da))),
  #           rawdata_df, colNames = FALSE, rowNames = FALSE)
}

for (i in 1:length(target_terms)) {
  
  saveWorkbook(get(paste0(target_terms[i],"_xlsx")),
               file = paste(number_of_polulation,"Gamma",paste("mr=",mr_gradient,collapse = " "),
                            target_terms[i],".xlsx",sep="_"))
}
pb$terminate()
