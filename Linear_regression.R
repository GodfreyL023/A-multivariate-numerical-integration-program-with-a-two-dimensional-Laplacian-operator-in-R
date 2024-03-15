remove(list=ls())
setwd("./Sample_data/")
library(ggplot2)
library(gridExtra)
library(xlsx)
library(openxlsx)
library(readxl)
library(paletteer)
library(tidyr)
library(lme4)
library(MASS)
library(nortest)
library(mgcv)
options(scipen = 100)
number_of_polulation <- 6
mr_gradient <- c(0,10^(-7:-2))
target_terms <- c("Biomass","Existence","Osci","Biomass_var",
                  "Existence_var","Osci_var")
target_terms2 <- c("Biomass","Species_diversity","Temporal_Instability","Patchiness",
                   "Species_aggregation","CV_of\ntemporal_instability")
par_df <- expand.grid(seq(0.25,2,by = 0.25),seq(0.25,2,by = 0.25))
colnames(par_df) <- c("beta1","beta2")
for (i in 1:length(target_terms)) {
  for (j in 1:length(mr_gradient)) {
    assign(paste("df",target_terms2[i], log10(mr_gradient[j]),sep = "_"),
           read_xlsx(paste(number_of_polulation,"Gamma",
                           paste("mr=",mr_gradient,collapse = " "),
                           target_terms[i],".xlsx",sep="_"), 
                     sheet = paste0("mr=",log10(mr_gradient[j])), col_names  = FALSE))
    assign(paste("df",target_terms2[i],log10(mr_gradient[j]), sep = "_"), 
           get(paste("df",target_terms2[i],log10(mr_gradient[j]), sep = "_"))[,-1])
    assign(paste("df",target_terms2[i],log10(mr_gradient[j]), sep = "_"), 
           cbind(par_df, get(paste("df",target_terms2[i],log10(mr_gradient[j]), sep = "_"))))
  }
}

mldf <- data.frame()
extended_par_df <- data.frame(rep(par_df[,1],each = 50), rep(par_df[,2],each = 50))
colnames(extended_par_df) <- c("beta1", "beta2")

for (i in 1:length(target_terms2)) {
  term_df <- data.frame()
  for (j in 1:length(mr_gradient)) {
    data <- get(paste("df",target_terms2[i],log10(mr_gradient[j]), sep = "_"))
    wide_data <- as.data.frame(t(data[,3:52]))
    colnames(wide_data) <- paste(data[,1],data[,2])
    long_data <- gather(wide_data, key = "pars", value = val)
    colnames(long_data)[2] <- target_terms2[i]
    long_data$Dr <- mr_gradient[j]
    if (j == 1){
      term_df <- cbind(extended_par_df, long_data[,2:3])
    } else {
      term_df <- rbind(term_df,cbind(extended_par_df, long_data[,2:3]))
    }
  }
  if (i == 1){
    mldf <- term_df
  } else {
    mldf[, ncol(mldf)+1] <- term_df[,3]
    colnames(mldf)[ncol(mldf)] <- target_terms2[i]
  }
}
# write.csv(mldf)

mldf_log <- mldf
for (i in 1:ncol(mldf_log)) {
  mldf_log[,i][mldf_log[,i] <= 0] <- min(mldf_log[,i][mldf_log[,i] > 0])/1000
  formula <- reformulate(c("poly(beta1, 1)",
                           "poly(beta2, 1)",
                           "poly(Dr, 1)",
                           "beta1 * Dr",
                           "beta2 * Dr",
                           "beta1 * beta2",
                           "beta1 * beta2 * Dr"), 
                         response = colnames(mldf_log)[i])
  bctrans <- boxcox(formula, data = mldf_log)
  lambda <- round(bctrans$x[which(bctrans$y == max(bctrans$y))],4)
  if (lambda == 0){
    mldf_log[,i] <- log(mldf_log[,i] + 1)
  } else {
    mldf_log[,i] <- (((mldf_log[,i] + 1)^lambda) - 1)/lambda
  }
}


formula <- reformulate(c("poly(beta1, 1)",
                         "poly(beta2, 1)",
                         "poly(Dr, 1)",
                         "beta1 * Dr",
                         "beta2 * Dr",
                         "beta1 * beta2",
                         "beta1 * beta2 * Dr"), response = target_terms2[3])
formula2 <- reformulate(c("beta1", 
                         "beta2", 
                         "Dr", 
                         "te(beta1, beta2)",
                         "te(beta1, Dr)",
                         "te(beta2, Dr)",
                         "te(beta1, beta2, Dr)"), response = target_terms2[1])
results <- lm(formula, data = mldf_log)
results2 <- gam(formula2, data = mldf_log)
hist(scale(results$residuals))
qqnorm(scale(results$residuals))
qqline(scale(results$residuals))
ks.test(scale(results$residuals), "pnorm")
hist(scale(results2$residuals))
qqnorm(scale(results2$residuals))
qqline(scale(results2$residuals))
ks.test(scale(results2$residuals), "pnorm")

summary(results)
summary(results2)


options(scipen = 0)
ggplot(mldf,aes(beta2, `Species aggregation`,group = Dr, color = Dr)) + 
  geom_point() +
  geom_smooth( se = FALSE) +
  scale_x_log10()+
  scale_y_log10()+
  ylim(0.1,1)+
  facet_wrap(~beta1)

