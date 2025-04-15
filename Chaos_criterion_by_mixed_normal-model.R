rm(list=ls())
path = "F:/AcademyData/Computer/rcode/ecology_model/GNM3_laplasian/local_dynamics/dynamic_type/dynamic_type_test2/"
setwd(path)
library(readxl)
library(mixtools)
library(nortest)
library(tidyr)
library(ggplot2)
library(paletteer)
mybeta <- seq(0.25,2,by=0.25)
my_sp_num <- 2:8
mydata <- data.frame()
for (i in 1:length(mybeta)) {
  for(j in 1:length(my_sp_num)){
    assign("data", read.csv(paste(paste("./GNM_beta", 
                                          round(mybeta[i],2), round(mybeta[i],2), 
                                          paste0("n=", my_sp_num[j]), sep = "_"), 
                                    "/", paste("GNM_beta", 
                                               round(mybeta[i],2), round(mybeta[i],2), 
                                               paste0("n=", my_sp_num[j]), "data.csv",sep = "_"), 
                                    sep = "")))
    data$beta <- mybeta[i]
    data$sp_num <- my_sp_num[j]
    mydata <- rbind(mydata, data)
  }
}

mixmodel <- normalmixEM(mydata$max_FTLE[which(mydata$dy_type==0&mydata$max_FTLE>0.002)], k=2)
normdf <- data.frame(Limit_cycles = dnorm(seq(min(mydata$max_FTLE[which(mydata$dy_type==0)]), 
                                          max(mydata$max_FTLE[which(mydata$dy_type==0)]), 
                                          by = 1e-4), 
                                      mean = mixmodel$mu[1], 
                                      sd = mixmodel$sigma[1]),
                     Chaos = dnorm(seq(min(mydata$max_FTLE[which(mydata$dy_type==0)]), 
                                          max(mydata$max_FTLE[which(mydata$dy_type==0)]), 
                                          by = 1e-4), 
                                      mean = mixmodel$mu[2], 
                                      sd = mixmodel$sigma[2]))
normdf <- gather(normdf, key = "distribution", value = "density")
normdf$variable <- rep(seq(min(mydata$max_FTLE[which(mydata$dy_type==0)]), 
                           max(mydata$max_FTLE[which(mydata$dy_type==0)]), 
                           by = 1e-4), 2)

# # 计算残差
# fitted_values <- apply(mixmodel$posterior, 1, function(posterior_row) {
#   sum(posterior_row * mixmodel$mu)
# })
# residuals <- mydata$max_FTLE[which(mydata$dy_type==0)] - fitted_values
# 
# # 绘制残差图
# plot(density(residuals))
# ks.test(residuals, "pnorm", mean=mean(residuals), sd=sd(residuals))
# ad.test(residuals)
# qqnorm(residuals)
# qqline(residuals, col="red")

calculate_posterior <- function(x, model) {
  lambda <- model$lambda
  mu <- model$mu
  sigma <- model$sigma
  
  # Calculate likelihood #
  likelihoods <- sapply(1:length(lambda), function(k) {
    dnorm(x, mean=mu[k], sd=sigma[k])
  })
  
  # Calculate posterior probability #
  posteriors <- likelihoods * lambda
  posteriors <- posteriors / sum(posteriors)
  
  return(posteriors)
}

mydata$dy_type[which(mydata$dy_type==0)] <- apply(calculate_posterior(mydata$max_FTLE[which(mydata$dy_type==0)], mixmodel), 
      MARGIN = 1, FUN = which.max) + 2

type_probs <- expand.grid(mybeta, my_sp_num)
type_probs[,3:6] <- 0
colnames(type_probs) <- c("beta", "sp_num", "1single_eq", "2multi_eq", "3limit_cyc", "4chaos")

find_criterion <- function(x){
  lambda <- mixmodel$lambda
  mu <- mixmodel$mu
  sigma <- mixmodel$sigma
  w_L1 <- dnorm(x, mu[1], sigma[1]) * lambda[1]
  w_L2 <- dnorm(x, mu[2], sigma[2]) * lambda[2]
  return((w_L1 - w_L2)/(w_L1 + w_L2))
}

FTLE_fig <- ggplot(mydata[which(mydata$dy_type >=3),], aes(x = max_FTLE)) +
  geom_histogram(binwidth = 0.002, fill = "grey", color = "black", alpha = 0.7) +
  geom_line(data = normdf, aes(x = variable, y = density*1.7,
                               group = distribution, color = distribution),
            linewidth = 1) +
  geom_vline(xintercept = uniroot(find_criterion, c(0.03, 0.07))$root, linetype = "dashed") +
  scale_y_continuous(name = "Frequency", sec.axis = sec_axis(~.*1.7, name="Density")) +
  scale_color_manual(values = c("yellow", "#2D3184FF"), labels = c("Chaos", "Limit Cycles")) +
  theme_classic() +
  labs(x = "Max FTLE")
ggsave("FTLE_distribution.eps", FTLE_fig, device = "eps", width = 5, height = 2.5)
ggsave("FTLE_distribution.png", FTLE_fig, width = 5, height = 2.5)

for(i in 1:length(mybeta)){
  for(j in 1:length(my_sp_num)){
    row_num <- which(mydata$beta==mybeta[i]&mydata$sp_num==my_sp_num[j])
    data <- mydata[row_num, ]
    myprobtable <- data.frame(table(data$dy_type)/sum(table(data$dy_type)))
    for (ii in 1:nrow(myprobtable)) {
      type_probs[which(type_probs$beta==mybeta[i]&type_probs$sp_num==my_sp_num[j]), (as.numeric(myprobtable[ii,1])+2)] <- myprobtable[ii,2]
    }
  }
}

custom_labeller <- function(variable, value) {
  if (variable == "d_type") {
    # 这里添加自定义的小标题，可以根据需要进行修改
    custom_titles <- c(
      "2multi_eq" = "Multi equilibrium",
      "3limit_cyc" = "Limit cycles",
      "4chaos" = "Chaos"
      # 添加其他类型和对应的小标题
    )
    return(custom_titles[value])
  }
  return(as.character(value))
}

type_prob_fig_df <- gather(type_probs[,4:6], key = "d_type", value = "prob")
type_prob_fig_df[,3] <- rep(mybeta, length(my_sp_num))
type_prob_fig_df[,4] <- rep(my_sp_num, each = length(mybeta))
colnames(type_prob_fig_df)[3:4] <- c("beta", "sp_num")
type_prob_fig <- ggplot(type_prob_fig_df, aes(beta, sp_num, fill = prob))+
  facet_wrap(~d_type, labeller = custom_labeller)+
  geom_tile(color = "black", size = 0.25) +
  scale_fill_gradient(low = "white", high = "orange", name = "", na.value = "white") +
  labs(x = expression(beta), y = "Species number") +
  theme_classic()
ggsave("Dynamic_type_prob.png", type_prob_fig, width = 8*0.8, height = 3*0.8)
ggsave("Dynamic_type_prob.eps", type_prob_fig, width = 8*0.8, height = 3*0.8)
