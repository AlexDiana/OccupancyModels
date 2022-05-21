library(here); library(MASS)
library(Rcpp); library(RcppArmadillo)

sourceCpp(here("MCMC/Functions","codecpp.cpp"))
source(here("Data/Simulated","functions.R"))
r_functions <- list.files(here("MCMC/Functions"))
r_functions <- r_functions[grep(".r", r_functions)]
sapply(r_functions, function(x) {
  source(here("MCMC/Functions",x))
})

# simulate data -----------------------------------------------------------

Y <- 5 # years
S_years <- rep(10000, Y)#floor(seq(2000, 2000, length.out = Y)) 
V_lambda <- 2 # number of visits

X <- expand.grid(seq(0, 1, length.out = sqrt(S_years[1])),
                 seq(0, 1, length.out = sqrt(S_years[1])))
# X <- cbind(runif(max(S_years), 0, 1), runif(max(S_years), 0, 1))
X <- scale(X)

# ggplot(data = NULL, aes(x = X[,1], y = X[,2])) + geom_point(size = .1)

mu_psi_true <- 0
beta_psi_tsp_true <- c(0, 0)
beta_psi_true <- c(1) # rep(0, 0)
sigma_gp <- .2
l_gp <- 2

K_l <- K(1:Y, 1:Y, sigma_gp, l_gp)
b_t_true <- mvrnorm(1, rep(0, Y), K_l)
# b_t_true <- mvrnorm(1, mu = rep(0, Y), Sigma = diag(0.2, nrow = Y))
#seq(-.5,.5,length.out = Y)#c(rep(-1,Y/2),rep(1,Y/2))
#mvrnorm(1, mu = rep(0, Y), Sigma = diag(0.1, nrow = Y))#seq(-.5,5,length.out = Y)

mu_p_true <- -1
beta_p_true <- rep(0, 0)#c(.5)##

sigma_s_true <- .5
l_s_true <- .2

sigma_eps_true <- .1

Sigma_X <- K2(X, X, sigma_s_true^2, l_s_true)
a_s_site_true <- rcpp_rmvnorm(1, Sigma_X, rep(0, nrow(X)))#rnorm(nrow(X), 0, sd = sigma_s_true)##########rep(0, nrow(X))###

# a_s_site_true <- c(rep(.5, max(S_years) * 75 / 100), rep(-.5, max(S_years) * 25 / 100))

# a_s_site_true <- a_s_site_true + rnorm(nrow(X), 0, sd = sigma_eps_true)

# a_s_site_true <- rnorm(nrow(X), 0, sd = sigma_eps_true)

a_s_site_true <- a_s_site_true - mean(a_s_site_true)
a_s_site_true <- a_s_site_true / 4

# plot effects
{
  ggplot() +
    geom_point(data = NULL, aes(x = X[,1],
                                y = X[,2], alpha = a_s_site_true), size = 3, shape = 15) +
    xlab("X") + ylab("Y") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          # legend.position = "none",
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          # panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
}

simulatedData <- simulateData(Y, S_years, V_lambda,
                              mu_psi_true, b_t_true, 
                              beta_psi_tsp_true, beta_psi_true,
                              X, a_s_site_true,
                              mu_p_true, beta_p_true)

colnames(simulatedData)[3:4] <- c("X1","X2")
if(length(beta_psi_true) > 0){
  colnames(simulatedData)[4 + seq_along(beta_psi_true)] <- paste0("OccCov",seq_along(beta_psi_true))  
}
if(length(beta_p_true) > 0){
  colnames(simulatedData)[4 + length(beta_psi_true) + seq_along(beta_p_true)] <- 
    paste0("CaptureCovariate",seq_along(beta_p_true))  
}

setwd(here("Data/Simulated"))

write.csv(simulatedData, file = "simdata.csv", row.names = F)
save(mu_psi_true, beta_psi_true, b_t_true, mu_p_true, a_s_site_true,
     beta_p_true, X, file = "simdata.rda")

