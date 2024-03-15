### This file exhibits a simple example of solving a multivariate differential system using deSolve

rm(list=ls())
path = "Set your own work direction"
dir.create(path)
setwd(path)
library(deSolve)
library(ggplot2)
library(patchwork)
library(tidyr)
library(paletteer)
library(xlsx)
library(openxlsx)

GNMmod <- function(Time, n, Pars) {
  with( as.list(Pars), {
    n[n < 0] <- 0
    dn <- r*n*(1 - n + bij%*%n - aij%*%(n^2))
    return(list(dn))
  })
}
model_type = "GNM"
number_of_polulation = 6
U <- c(2,2)
iterations = 100
nn <- number_of_polulation^2
times <- seq(0, 100, by = 0.01)
aij <- matrix(rgamma(nn,shape = 1,scale = U[1]),
              nrow = number_of_polulation,ncol = number_of_polulation)
bij <- matrix(rgamma(nn,shape = 1,scale = U[2]),
              nrow = number_of_polulation,ncol = number_of_polulation)
diag(aij) <- 0
diag(bij) <- 0
r = runif(number_of_polulation,0.5,1)
pars <- c(aij,bij,r)
NN <- runif(number_of_polulation,0,1)
out <- ode(NN, times, GNMmod, pars)
out_data <- gather(data.frame(out[,2:ncol(out)]),key = "species", value = "abundance")
out_data$time <- times
local_dynamic <- ggplot(out_data,aes(time,abundance,group = species))+
  geom_line(aes(color = species),show.legend = FALSE,linewidth = 0.5) +
  labs(x = "Time", y = "Abundance") +
  scale_color_manual(values = paletteer_c("grDevices::Blue-Yellow", number_of_polulation+2)) +
  theme_classic()
local_dynamic

type = "Chaos4"
if (grepl("Stable", type)){
  local_dynamic <- local_dynamic + ggtitle("e") + 
    theme(plot.title = element_text(size = 18, face = "bold"))
} else if(grepl("Osci", type)){
  local_dynamic <- local_dynamic + ggtitle("f") + 
    theme(plot.title = element_text(size = 18, face = "bold"))
} else{
  local_dynamic <- local_dynamic + ggtitle("g") + 
    theme(plot.title = element_text(size = 18, face = "bold"))
}
ggsave(paste0("./new3/", type, ".pdf"), 
       local_dynamic, width = 3, height = 2.75)
mybook <- createWorkbook()
my_aij <- addWorksheet(mybook, sheetName = "aij")
writeData(mybook, my_aij, aij, colNames = FALSE, rowNames = FALSE)
my_bij <- addWorksheet(mybook, sheetName = "bij")
writeData(mybook, my_bij, bij, colNames = FALSE, rowNames = FALSE)
my_r <- addWorksheet(mybook, sheetName = "r")
writeData(mybook, my_r, r, colNames = FALSE, rowNames = FALSE)
my_out <- addWorksheet(mybook, sheetName = "out")
writeData(mybook, my_out, out_data, colNames = TRUE, rowNames = FALSE)
saveWorkbook(mybook, file = paste0("./new3/", type, ".xlsx"))

