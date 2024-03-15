### In many special models, the space always change over time under the combined action of instable dynamics of the model and diffusion.
### The change of space through time can be visualized using this file.

library(ggplot2)
library(gridExtra)
library(ggpubr)
rm(list=ls())
setwd("F:/AcademyData/Computer/rcode/ecology_model/GNM3_laplasian/preliminary_experiment/case46/")
options(scipen = 100)
library(readxl)
iterations = 100
N = 50
number_of_population <- 6
U <- c(2,2)
dimension_of_space = N
dd = dimension_of_space^2
Da = 1e-3
Pdata <- read_excel(paste(N, number_of_population, "Gamma", U[1], U[2],
                          "mr=", Da, ".csv", sep ="_"), sheet = "Population Data")
#Cdata <- read_excel(paste("GNM_20_0_5_mr=",mr,".csv",sep="_"), sheet = "Community Data")

myc <- seq(from = 86, by = 2, length.out = 8)

data <- data.frame(Pdata[,myc[1]],
                   1:dimension_of_space,
                   rep(1:dimension_of_space,each=dimension_of_space))
names(data) <- c("abundance","x","y")
data$abundance <- (data$abundance-min(data$abundance))/(max(data$abundance)-min(data$abundance))
p1 <- ggplot(data,aes(x,y,fill = abundance))+
  geom_tile(color = "white", size = 0.25) +
  scale_fill_gradient2(low = "#2D3184FF", mid = "#1296AEFF", high = "yellow",
                       midpoint = 0.4,limits = c(0,1), name = "Relative\nsize") +
  labs(x = "", y = "") +
  ggtitle(paste0("Time = ", myc[1])) +
  theme_classic()

data <- data.frame(Pdata[,myc[2]],
                   1:dimension_of_space,
                   rep(1:dimension_of_space,each=dimension_of_space))
names(data) <- c("abundance","x","y")
data$abundance <- (data$abundance-min(data$abundance))/(max(data$abundance)-min(data$abundance))
p2 <- ggplot(data,aes(x,y,fill = abundance))+
  geom_tile(color = "white", size = 0.25) +
  scale_fill_gradient2(low = "#2D3184FF", mid = "#1296AEFF", high = "yellow",
                       midpoint = 0.4,limits = c(0,1), name = "Relative\nsize") +  labs(x = "", y = "") +
  ggtitle(paste0("Time = ", myc[2])) +
  theme_classic()

data <- data.frame(Pdata[,myc[3]],
                   1:dimension_of_space,
                   rep(1:dimension_of_space,each=dimension_of_space))
names(data) <- c("abundance","x","y")
data$abundance <- (data$abundance-min(data$abundance))/(max(data$abundance)-min(data$abundance))
p3 <- ggplot(data,aes(x,y,fill = abundance))+
  geom_tile(color = "white", size = 0.25) +
  scale_fill_gradient2(low = "#2D3184FF", mid = "#1296AEFF", high = "yellow",
                       midpoint = 0.4,limits = c(0,1), name = "Relative\nsize") +  labs(x = "", y = "") +
  ggtitle(paste0("Time = ", myc[3])) +
  theme_classic()

data <- data.frame(Pdata[,myc[4]],
                   1:dimension_of_space,
                   rep(1:dimension_of_space,each=dimension_of_space))
names(data) <- c("abundance","x","y")
data$abundance <- (data$abundance-min(data$abundance))/(max(data$abundance)-min(data$abundance))
p4 <- ggplot(data,aes(x,y,fill = abundance))+
  geom_tile(color = "white", size = 0.25) +
  scale_fill_gradient2(low = "#2D3184FF", mid = "#1296AEFF", high = "yellow",
                       midpoint = 0.4,limits = c(0,1), name = "Relative\nsize") +  labs(x = "", y = "") +
  ggtitle(paste0("Time = ", myc[4])) +
  theme_classic()

data <- data.frame(Pdata[,myc[5]],
                   1:dimension_of_space,
                   rep(1:dimension_of_space,each=dimension_of_space))
names(data) <- c("abundance","x","y")
data$abundance <- (data$abundance-min(data$abundance))/(max(data$abundance)-min(data$abundance))
p5 <- ggplot(data,aes(x,y,fill = abundance))+
  geom_tile(color = "white", size = 0.25) +
  scale_fill_gradient2(low = "#2D3184FF", mid = "#1296AEFF", high = "yellow",
                       midpoint = 0.4,limits = c(0,1), name = "Relative\nsize") +  labs(x = "", y = "") +
  ggtitle(paste0("Time = ", myc[5])) +
  theme_classic()

data <- data.frame(Pdata[,myc[6]],
                   1:dimension_of_space,
                   rep(1:dimension_of_space,each=dimension_of_space))
names(data) <- c("abundance","x","y")
data$abundance <- (data$abundance-min(data$abundance))/(max(data$abundance)-min(data$abundance))
p6 <- ggplot(data,aes(x,y,fill = abundance))+
  geom_tile(color = "white", size = 0.25) +
  scale_fill_gradient2(low = "#2D3184FF", mid = "#1296AEFF", high = "yellow",
                       midpoint = 0.4,limits = c(0,1), name = "Relative\nsize") +  labs(x = "", y = "") +
  ggtitle(paste0("Time = ", myc[6])) +
  theme_classic()

data <- data.frame(Pdata[,myc[7]],
                   1:dimension_of_space,
                   rep(1:dimension_of_space,each=dimension_of_space))
names(data) <- c("abundance","x","y")
data$abundance <- (data$abundance-min(data$abundance))/(max(data$abundance)-min(data$abundance))
p7 <- ggplot(data,aes(x,y,fill = abundance))+
  geom_tile(color = "white", size = 0.25) +
  scale_fill_gradient2(low = "#2D3184FF", mid = "#1296AEFF", high = "yellow",
                       midpoint = 0.4,limits = c(0,1), name = "Relative\nsize") +  labs(x = "", y = "") +
  ggtitle(paste0("Time = ", myc[7])) +
  theme_classic()

data <- data.frame(Pdata[,myc[8]],
                   1:dimension_of_space,
                   rep(1:dimension_of_space,each=dimension_of_space))
names(data) <- c("abundance","x","y")
data$abundance <- (data$abundance-min(data$abundance))/(max(data$abundance)-min(data$abundance))
p8 <- ggplot(data,aes(x,y,fill = abundance))+
  geom_tile(color = "white", size = 0.25) +
  scale_fill_gradient2(low = "#2D3184FF", mid = "#1296AEFF", high = "yellow",
                       midpoint = 0.4,limits = c(0,1), name = "Relative\nsize") +  labs(x = "", y = "") +
  ggtitle(paste0("Time = ", myc[8])) +
  theme_classic()

p9 <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow = 4,ncol=2,common.legend = TRUE, legend="right")
p10 <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow=2,ncol=4,common.legend = TRUE, legend="right")

p9
ggsave(paste(N, number_of_population, "Gamma", U[1], U[2],Da,"_snapshot.pdf",sep="_"),p10, width = 8,height = 4)
ggsave(paste(N, number_of_population, "Gamma", U[1], U[2],Da,"_snapshot2.pdf",sep="_"),p9, width = 4.5,height = 8)


