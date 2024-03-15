### Use this program to redrawn the graphs using the data obtained in "Laplacian.R"
rm(list=ls())

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(readxl)
library(tidyr)
library(paletteer)
setwd("set your own work direction")
dir.create("./redrawn")
options(scipen = 100)
iterations = 100
N = 50
number_of_population <- 10
U <- c(2, 2)
dimension_of_space = N
dd = dimension_of_space^2
Da = c(0, 1e-6, 1e-5, 1e-4)
for (i in 1:length(Da)) {
  assign(paste0("Pdata", i), read_excel(paste(N, number_of_population, "Gamma", U[1], U[2],
                                              "mr=", Da[i], ".csv", sep ="_"), sheet = "Population Data"))
  assign(paste0("Cdata", i), read_excel(paste(N, number_of_population, "Gamma", U[1], U[2],
                                              "mr=", Da[i], ".csv", sep ="_"), sheet = "Community Data"))
}
sample_size <- ceiling(dd/200) 
sample_grids <- sample(1:dd, sample_size)

line_colors <- paletteer_c("grDevices::Blue-Yellow",sample_size) #"grDevices::Viridis"
line_colors <- line_colors[sample(1:50, 50)]
for (i in 1:length(Da)) {
  if(Da[i] < 0.01 & Da[i] != 0){
    myd <- format(Da[i],scientific = TRUE)
  }else{
    myd = Da[i]
  }
  Pdata <- get(paste0("Pdata",i))
  data <- data.frame(Pdata[,iterations+1],
                     1:dimension_of_space,
                     rep(1:dimension_of_space,each=dimension_of_space))
  names(data) <- c("abundance","x","y")
  data$abundance <- (data$abundance)/(max(data$abundance))
  p2 <- ggplot(data,aes(x,y,fill = abundance))+
    geom_tile(color = "white", size = 0.25) +
    scale_fill_gradient2(low = "#2D3184FF", mid = "#1296AEFF", high = "yellow",
                         midpoint = 0.4,limits = c(0,1), name = "Relative\nsize") +
    labs(x = "x", y = "y") +
    theme_classic() +
    theme(legend.key.size = unit(0.4, "cm"), 
          legend.text = element_text(size = 8), 
          legend.title = element_text(size = 10))
  assign(paste0("spaceP", i), p2)
  
  data <- gather(data.frame(t(Pdata[sample_grids, ])), key = "grid", value = "abundance")
  data$times <- rep(1:(iterations+1), sample_size)
  p1 <- ggplot(data,aes(times,abundance),group = grid)+
    geom_line(aes(color = grid), 
              show.legend = FALSE, 
              linewidth = 0.7) +
    labs(x = "Time", y = "Population size") +
    ggtitle(bquote(italic("D") == .(myd))) +
    scale_color_manual(values = line_colors) +
    theme_classic()
  assign(paste0("population", i), p1)
  
  Cdata <- get(paste0("Cdata",i))
  data <- data.frame(Cdata[,iterations+1],
                     1:dimension_of_space,
                     rep(1:dimension_of_space,each=dimension_of_space))
  names(data) <- c("abundance","x","y")
  data$abundance <- (data$abundance)/(max(data$abundance))
  p5 <- ggplot(data,aes(x,y,fill = abundance))+
    geom_tile(color = "white", size = 0.25) +
    scale_fill_gradient2(low = "#2D3184FF", mid = "#1296AEFF", high = "yellow",
                         midpoint = 0.4,limits = c(0,1), name = "Relative\nsize") +
    labs(x = "x", y = "y") +
    theme_classic() + 
    theme(legend.key.size = unit(0.4, "cm"), 
          legend.text = element_text(size = 8), 
          legend.title = element_text(size = 10))
  assign(paste0("spaceC", i), p5)
  
  data <- gather(data.frame(t(Cdata[sample_grids,])), key = "grid", value = "abundance")
  data$times <- rep(1:(iterations+1),sample_size)
  p4 <- ggplot(data,aes(times,abundance),group = grid)+
    geom_line(aes(color = grid),show.legend = FALSE,linewidth = 0.7) +
    labs(x = "Time", y = "Total population size") +
    ggtitle(bquote(italic("D") == .(myd))) +
    scale_color_manual(values = line_colors) +
    theme_classic()
  assign(paste0("Total_population", i), p4)
}

Pfig <- grid.arrange(population1, population2, population3, population4, 
                     spaceP1, spaceP2, spaceP3, spaceP4, nrow = 2, heights = c(1.2,1))
ggsave(paste("./redrawn/", number_of_population, "Gamma", U[1], U[2],
             "Population.pdf",sep ="_"),
       Pfig,dpi=1000,width = 11,height = 4.5)

Cfig <- grid.arrange(Total_population1, Total_population2, Total_population3, Total_population4, 
                     spaceC1, spaceC2, spaceC3, spaceC4, nrow = 2, heights = c(1.2,1))
ggsave(paste("./redrawn/", number_of_population, "Gamma", U[1], U[2],
             "Community.pdf",sep ="_"),
       Cfig,dpi=1000,width = 11,height = 4.5)
