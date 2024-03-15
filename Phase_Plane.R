### The dynamics of multivariate systems are not intuitive, and most non-linear systems are not analytically solvable.
### A useful way to begin to understand a system is to perform the phase plane analysis on the simple cases of the multi-dimension system.
### This file presents a simple program to do the phase plane analysis for the model in our article.

rm(list=ls())

setwd("F:/AcademyData/Computer/rcode/ecology_model/GNM3_laplasian/phase")
library(deSolve)
library(phaseR)

GNM <- function(t,y,parameters){
  x <- y[1]
  y <- y[2]
  a1 <- parameters[1]
  a2 <- parameters[2]
  b1 <- parameters[3]
  b2 <- parameters[4]
  dy <- numeric(2)
  dy[1] <- x*( 1 - x - a1*y^2 + b1*y)
  dy[2] <- y*( 1 - y - a2*x^2 + b2*x)
  list(dy)
}

# NN <- runif(2, 0, 1)
# Times <- seq(0,100,by = 0.1)
Pars <- matrix(c(c(0.5,0.5,0.5,0.5), 
                 c(1,1,1,1), 
                 c(1,1,1,1.5),
                 c(2,2,2,2)), nrow = 4)
# out <- ode(NN, Times, GNM, Pars)
# matplot(out[,2:3], type = "l")
# for (i in 1:ncol(Pars)) {
#   
# }
figure_no <- 4
alphabet <- c("a", "b", "c", "d")
pp <- pdf(paste0("a1=",Pars[1,figure_no],"_a2=",Pars[2,figure_no],"_b1=",Pars[3,figure_no],
                 "_b2=",Pars[4,figure_no],"phase plainnn.pdf"),height = 3*1.25,width = 3*1.25)

flowField(GNM,xlim = c(0,2),ylim = c(0,2),parameters = Pars[,figure_no],
          points = 20,add = FALSE,col = "#808080")
grid()

nullclines(GNM,xlim = c(0,2),ylim = c(0,2),parameters = Pars[,figure_no],
           col=c("#2D3184FF","yellow"),lwd=2,add.legend = FALSE)
legend("topright",c("y","x"),lty = c(1,1),lwd=c(2,2),col=c("#2D3184FF","yellow"))

##
Interval_step <- seq(0,2,by = 0.1)
Interval_step[1] <- 0.001
ff <- function(x,parameters = Pars[,figure_no]){
  a1 <- parameters[1]
  a2 <- parameters[2]
  b1 <- parameters[3]
  b2 <- parameters[4]
  y <- 1 - a2*x^2 + b2*x
  z <- 1 - x - a1*y^2 + b1*y
  return(z)
}

for (i in 1:(length(Interval_step)-1)) {
  tryCatch({
    xx <- uniroot(ff,interval = c(Interval_step[i],Interval_step[i+1]))$root
    yy <- 1 - Pars[2, figure_no]*xx^2 + Pars[4, figure_no]*xx
    tp <- stability(GNM,ystar = c(xx,yy),parameters = Pars[,figure_no])$classification
    if(tp == "Stable node" |tp == "Stable focus"){
      points(xx,yy,pch = 19)
    }
  },error = function(e){})
}
mtext(alphabet[figure_no], adj = 0, cex = 1.5, font = 2)
dev.off()

