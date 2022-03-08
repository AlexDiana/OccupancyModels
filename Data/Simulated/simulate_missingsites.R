
# simulate data missing sites -----------------------------------------------------------

Y <- 40 # years
S <- 10000 # sites
V_lambda <- .5 # number of visits
p_sitespresence <- .05

X <- expand.grid(seq(0, 1, length.out = sqrt(S)),
                 seq(0, 1, length.out = sqrt(S)))
# X <- cbind(runif(S, 0, 1), runif(S, 0, 1))
X <- scale(X)

mu_psi_true <- 0
beta_psi_tsp_true <- c(0, 0)
beta_psi_true <- c()
sigma_gp <- .2
l_gp <- 1

K_l <- K(1:Y, 1:Y, sigma_gp, l_gp)
b_t_true <- mvrnorm(1, rep(0, Y), K_l)

mu_p_true <- -1
beta_p_true <- rep(1, 0)

sigma_s_true <- .5
l_s_true <- .3

sigma_eps_true <- .5

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
                                y = X[,2], alpha = a_s_site_true), size = 1.5, shape = 15) +
    xlab("X") + ylab("Y") + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          # legend.position = "none",
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          # panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
}

# a_s_site_true <- rep(a_s_site_true_base, Y)

simulatedData <- simulateDataMissingSites(Y, S, V_lambda,
                                          mu_psi_true, b_t_true, 
                                          beta_psi_tsp_true, beta_psi_true,
                                          X, a_s_site_true,
                                          mu_p_true, beta_p_true, p_sitespresence)


colnames(simulatedData)[3:4] <- c("X1","X2")
if(length(beta_psi_true) > 0){
  colnames(simulatedData)[4 + seq_along(beta_psi_true)] <- paste0("OccCov",seq_along(beta_psi_true))  
}
if(length(beta_p_true) > 0){
  colnames(simulatedData)[4 + length(beta_psi_true) + seq_along(beta_p_true)] <- 
    paste0("CaptureCovariate",seq_along(beta_p_true))  
}

setwd("/cluster/home/osr/ad625/FastOccupancy/Simulated Datasets")

# write.csv(simulatedData, file = paste0("simdata_",S_years[1],"sites.csv"), row.names = F)
# save(mu_psi_true, beta_psi_true, b_t_true, mu_p_true, a_s_site_true,
# beta_p_true, X, file = paste0("simdata_",S_years[1],"sites.rda"))

write.csv(simulatedData, file = "simdata_spatial3.csv", row.names = F)
save(mu_psi_true, beta_psi_true, b_t_true, mu_p_true, a_s_site_true,
     beta_p_true, X, file = "simdata_spatial3.rda")

length(unique(simulatedData$Site))
