library(MASS); library(dplyr); 
library(coda); library(FastGP); library(here)
library(reshape2); library(ggplot2); 
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
S_years <- rep(3600, Y)#floor(seq(2000, 2000, length.out = Y)) 
V_lambda <- 2 # number of visits

X <- expand.grid(seq(0, 1, length.out = sqrt(S_years[1])),
                 seq(0, 1, length.out = sqrt(S_years[1])))
X <- scale(X)



# ggplot(data = NULL, aes(x = X[,1], y = X[,2])) + geom_point(size = .1)

mu_psi_true <- 0
beta_psi_tsp_true <- c(0, 0)
beta_psi_true <- rep(0, 0)
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
l_s_true <- .25

sigma_eps_true <- .1

Sigma_X <- K2(X, X, sigma_s_true^2, l_s_true)
a_s_site_true <- rcpp_rmvnorm(1, Sigma_X, rep(0, nrow(X)))#rnorm(nrow(X), 0, sd = sigma_s_true)##########rep(0, nrow(X))###

a_s_site_true <- a_s_site_true - mean(a_s_site_true)

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

# idx_training <- seq(2, nrow(X), 2)
gridLength <- sqrt(nrow(X))
idx_training <- rep(NA, gridLength * gridLength / 2)
l <- 0
for (i in 1:gridLength) {
  if(i %% 2 == 0){ # even
    idx_training[(i - 1) * (gridLength / 2) + 1:(gridLength / 2)] <- (i - 1) * gridLength + seq(2, gridLength, 2)
  } else {
    idx_training[(i - 1) * (gridLength / 2) + 1:(gridLength / 2)] <- (i - 1) * gridLength + seq(2, gridLength, 2) - 1
  }
}

X_training <- X[idx_training,]
X_test <- X[-idx_training,]

idx_color <- rep(0, nrow(X))
idx_color[idx_training] <- 1

ggplot(data = NULL, aes(x = X[,1], y = X[,2],
                        color = idx_color)) + geom_point(size = .1)

a_s_training <- a_s_site_true[idx_training]
a_s_test <- a_s_site_true[-idx_training]

S_years <- rep(S_years[1] * .5, Y)

simulatedData <- simulateData(Y, S_years, V_lambda,
                              mu_psi_true, b_t_true, 
                              beta_psi_tsp_true, beta_psi_true,
                              X_training, a_s_training,
                              mu_p_true, beta_p_true)

colnames(simulatedData)[3:4] <- c("X1","X2")

setwd(here("Data/Simulated/GPcomparison"))

write.csv(simulatedData, file = "simdata.csv", row.names = F)
save(mu_psi_true, beta_psi_true, b_t_true, mu_p_true, a_s_training, a_s_test,
     beta_p_true, X_training, X_test, file = "simdata.rda")



# run models ---------

realData <- F

# input
{
  data_file <- here("Data/Simulated/GPcomparison","simdata.csv")
  
  data <- read.csv(file = data_file, stringsAsFactors = F)
  index_year <- 1
  index_site <- 2
  index_occ <- 5
  index_spatial_x <- 3
  index_spatial_y <- 4
  covariates_psi_text <- "0"#  "7"
  covariates_p_text <- "0"  #"5-6"
  
}

# prior parameters
{
  prior_psi <- .1
  sigma_psi <- 2
  
  prior_p <- .1
  sigma_p <- 2
  
  # buildSpatialGrid(data$X1, data$X2, .1)
}

# mcmc parameter
{
  nchain <- 1
  nburn <- 2000#15000
  niter <- 2000#20000
}

# SOR ------

usingSpatial <- T
usingYearDetProb <- F
spatialApprox <- "SoR"
storeRE <- F

grid_values <- c(.6, .4, .3, .2)
MSE_grid <- rep(NA, length(grid_values))

for (idx_grid in seq_along(grid_values)) {
  
  gridStep <- grid_values[idx_grid]
  
  maxPoints <- 5
  
  modelResults_MCMC <- runModel(data, 
                                index_year, # index of the column with year
                                index_site, # index of the column with site
                                index_occ, # index of the column with captures
                                index_spatial_x, # index of the column with the x coordinate
                                index_spatial_y, # index of the column with the y coordinate
                                covariates_psi_text, # index of numerical covariates for psi
                                covariates_p_text, # index of numerical covariates for p
                                prior_psi, # prior mean of occupancy probability
                                sigma_psi, # just as in the notes
                                prior_p, 
                                sigma_p, 
                                usingYearDetProb,
                                usingSpatial,
                                spatialApprox,
                                maxPoints,
                                gridStep, # width of the spatial grid
                                storeRE = F,
                                nchain, # number of chains
                                nburn, # burn-in iterations
                                niter) # number of iterations
  
  # compare predictions -----
  
  X_tilde <- modelResults_MCMC$dataCharacteristics$X_tilde
  X_centers <- nrow(X_tilde)
  years <- modelResults_MCMC$dataCharacteristics$Years
  Y <- length(years)
  
  beta_psi_output <- modelResults_MCMC$modelOutput$beta_psi_output
  beta_psi_output <- apply(beta_psi_output, 3, c)
  
  l_s_output <- modelResults_MCMC$modelOutput$l_s_output
  l_s_mean <- mean(l_s_output)
  
  a_s_pred_test_output <- sapply(1:niter, function(i){
    print(i)
    K_staru <- K2(X_test, X_tilde, 1, l_s_output[i])
    inv_K_uu <- solve(K2(X_tilde, X_tilde, 1, l_s_output[i]) + 
                        diag(exp(-10), nrow = nrow(X_tilde)))
    X_coeff <- K_staru %*% inv_K_uu
    
    X_coeff %*% beta_psi_output[i, Y + seq_len(X_centers)]
  })
  
  # K_staru <- K2(X_test, X_tilde, 1, l_s_mean)
  # inv_K_uu <- solve(K2(X_tilde, X_tilde, 1, l_s_mean) +
  #                     diag(exp(-10), nrow = nrow(X_tilde)))
  # X_coeff <- K_staru %*% inv_K_uu
  # 
  # beta_mean <- apply(beta_psi_output, 2, mean)
  # 
  # a_s_pred_test <- X_coeff %*% beta_mean[Y + seq_len(X_centers)]
  
  a_s_pred_test <- apply(a_s_pred_test_output, 1, mean)
  
  MSE_grid[idx_grid] <-  mean(abs(a_s_pred_test - a_s_test))
  
  # qplot(a_s_pred_test - a_s_test)
  # tracePlot_OccupancySiteEffect(modelResults_MCMC, 35)
  
}

qplot(grid_values, MSE_grid)

# SOD ------------

usingSpatial <- T
usingYearDetProb <- F
spatialApprox <- "SoD"
storeRE <- F

grid_values <- c(.6, .4, .3, .2)
MSE_grid <- rep(NA, length(grid_values))

for (idx_grid in seq_along(grid_values)) {
  
  gridStep <- grid_values[idx_grid]
  
  maxPoints <- NA
  
  modelResults_MCMC <- runModel(data, 
                                index_year, # index of the column with year
                                index_site, # index of the column with site
                                index_occ, # index of the column with captures
                                index_spatial_x, # index of the column with the x coordinate
                                index_spatial_y, # index of the column with the y coordinate
                                covariates_psi_text, # index of numerical covariates for psi
                                covariates_p_text, # index of numerical covariates for p
                                prior_psi, # prior mean of occupancy probability
                                sigma_psi, # just as in the notes
                                prior_p, 
                                sigma_p, 
                                usingYearDetProb,
                                usingSpatial,
                                spatialApprox,
                                maxPoints,
                                gridStep, # width of the spatial grid
                                storeRE = F,
                                nchain, # number of chains
                                nburn, # burn-in iterations
                                niter) # number of iterations
  
  # compare predictions -----
  
  X_tilde <- modelResults_MCMC$dataCharacteristics$X_tilde
  X_centers <- nrow(X_tilde)
  years <- modelResults_MCMC$dataCharacteristics$Years
  Y <- length(years)
  
  beta_psi_output <- modelResults_MCMC$modelOutput$beta_psi_output
  beta_psi_output <- apply(beta_psi_output, 3, c)
  
  indexesSite <- findClosestPoint(as.matrix(X_test), X_tilde)
  
  a_s_pred_test <- sapply(1:nrow(X_test), function(i){
    mean(beta_psi_output[, Y + indexesSite[i]]
    )
  })
  
  MSE_grid[idx_grid] <-  mean(abs(a_s_pred_test - a_s_test))
  
}
