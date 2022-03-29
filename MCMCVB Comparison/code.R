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

S_grid <- c(500, 1000, 2500, 5000)

beta_effect_MCMC_all <- matrix(NA, 3, length(S_grid))
beta_effect_VB_all <- matrix(NA, 3, length(S_grid))

for (idx_s in seq_along(S_grid)) {
  
  S_current <- S_grid[idx_s]
  
  # simulate data
  {
    Y <- 10 # years
    S_years <- rep(S_current, Y)#floor(seq(2000, 2000, length.out = Y)) 
    V_lambda <- 2 # number of visits
    
    X <- expand.grid(seq(0, 1, length.out = sqrt(S_years[1])),
                     seq(0, 1, length.out = sqrt(S_years[1])))
    X <- scale(X)
    
    mu_psi_true <- 0
    beta_psi_tsp_true <- c(0, 0)
    beta_psi_true <- c(1) # rep(0, 0)
    sigma_gp <- .2
    l_gp <- 2
    
    K_l <- K(1:Y, 1:Y, sigma_gp, l_gp)
    b_t_true <- mvrnorm(1, rep(0, Y), K_l)
    
    mu_p_true <- -1
    beta_p_true <- rep(0, 0)
    
    sigma_s_true <- .5
    l_s_true <- .25
    
    sigma_eps_true <- .1
    
    a_s_site_true <- rnorm(nrow(X), 0, sd = sigma_eps_true)
    
    a_s_site_true <- a_s_site_true - mean(a_s_site_true)
    
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
    
    setwd(here("Data/Simulated/MCMCVB"))
    
    data_name <- paste0("simdata_",S_current)
    
    write.csv(simulatedData, file = paste0(data_name,".csv"), row.names = F)
    save(mu_psi_true, beta_psi_true, b_t_true, mu_p_true, a_s_site_true,
         beta_p_true, X, file = paste0(data_name,".rda"))
  }
  
  # data input
  {
    # input
    {
      
      data_file <- here("Data/Simulated/MCMCVB",paste0(data_name,".csv"))
      
      data <- read.csv(file = data_file, stringsAsFactors = F)
      index_year <- 1
      index_site <- 2
      index_occ <- 6
      index_spatial_x <- 3
      index_spatial_y <- 4
      covariates_psi_text <- "5"#  "7"
      covariates_p_text <- "0"  #"5-6"
      
    }
    
    # prior parameters
    {
      prior_psi <- .1
      sigma_psi <- 2
      
      prior_p <- .1
      sigma_p <- 2
      
      gridStep <-  .2025
    }
    
    
    usingSpatial <- F
    usingYearDetProb <- F
    spatialApprox <- ""
    maxPoints <- 5
    storeRE <- F
  }
  
  # mcmc
  {
    sourceCpp(here("MCMC/Functions","codecpp.cpp"))
    r_functions <- list.files(here("MCMC/Functions"))
    r_functions <- r_functions[grep(".r", r_functions)]
    sapply(r_functions, function(x) {
      source(here("MCMC/Functions",x))
    })
    
    # mcmc parameter
    {
      nchain <- 1
      nburn <- 2500
      niter <- 2500
    }
    
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
  }
  
  # vb
  {
    sourceCpp(here("VB/Functions","codecpp.cpp"))
    r_functions <- list.files(here("VB/Functions"))
    r_functions <- r_functions[grep(".r", r_functions)]
    sapply(r_functions, function(x) {
      source(here("VB/Functions",x))
    })
    
    # algorithm parameter
    {
      nruns <- 1
      tol <- .01#15000
      maxit <- 1000
      
      verbose <- T
      
    }
    
    modelResults_VB <- runModel(data, 
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
                                gridStep, # width of the spatial grid
                                storeRE = F,
                                nruns, # number of chains
                                tol, # burn-in iterations
                                maxit) # number of iterations
  }
  
  # save results
  {
    
    setwd("~/OccupancyModels/MCMCVB Comparison")
    
    # years ---
    {
      
      years <- modelResults_VB$dataCharacteristics$Years
      Y <- length(years)
      
      # mcmc 
      {
        X_tilde <- modelResults_MCMC$dataCharacteristics$X_tilde
        if(usingSpatial){
          X_centers <- nrow(X_tilde)  
        } else {
          X_centers <- 0
        }
        
        beta_psi_output <- modelResults_MCMC$modelOutput$beta_psi_output
        beta_psi_output <- apply(beta_psi_output, 3, c)
        if(!usingSpatial){
          # a_s_unique_output <- modelResults$a_s_unique_output
          # a_s_unique_output <- apply(a_s_unique_output, 3, c)
        }
        
        CI_yearseffect_MCMC <- sapply(1:Y, function(j) {
          
          if(usingSpatial){
            yearEffect <- logit(
              beta_psi_output[,j] #+
              # apply(modelResults$a_s_unique_output[1,,], 1, mean) +
              # apply(beta_psi_output[,Y + 1:X_centers], 1, mean)
            )
          } else {
            yearEffect <- logit(
              # beta_psi_output[, 1] +
              beta_psi_output[, j] 
              # X_psi_yearcov[j] * beta_psi_output[,1 + Y + X_centers + 1] + 
              # apply(a_s_unique_output, 1, mean)
            )
          }
          
          c(quantile(yearEffect, probs = c(0.025, 0.975)), mean(yearEffect))
          
        })
      }
      
      # vb
      {
        mu_beta_psi <- modelResults_VB$modelOutput$mu_beta_psi
        Sigma_beta_psi <- modelResults_VB$modelOutput$Sigma_beta_psi
        
        CI_yearseffect_VB <- sapply(1:Y, function(j){
          logit(mu_beta_psi[j] + c(-1,1,0) * 1.96 * sqrt(Sigma_beta_psi[j,j]))
        })
        
      }
      
      prob_true <- logit(mu_psi_true + b_t_true #+ 
                         # beta_psi_t_true * sort(unique(modelResults$dataCharacteristics$X_psi_yearcov)) +
                         # mean(a_s_site_true)
      )
      
      y_plot <- ggplot(data = NULL, aes(x = years,
                                        y = CI_yearseffect_MCMC[3,],
                                        ymin = CI_yearseffect_MCMC[1,],
                                        ymax = CI_yearseffect_MCMC[2,])) + 
        # geom_point() + 
        geom_errorbar(alpha = .4, size = 1) + 
        geom_line(data = NULL, aes(x = years, y = prob_true), color = "red")  +
        # geom_point(data = NULL, aes(x = years,
        # y = CI_yearseffect_VB[3,]), color = "black") +
        geom_errorbar(data = NULL, aes(x = years,
                                       ymin = CI_yearseffect_VB[1,],
                                       ymax = CI_yearseffect_VB[2,]), color = "black", size = 1) +
        geom_point(data = NULL, aes(x = years, y = prob_true), size = 2, color = "red") +
        ylab("Occupancy probability") + scale_x_continuous(name = "Year", breaks = years) +
        theme(plot.title = element_text(hjust = 0.5, size = 20),
              axis.title = element_text(size = 20, face = "bold"),
              axis.text = element_text(size = 13, face = "bold", angle = 90),
              panel.grid.major = element_line(colour="grey", size=0.015),
              panel.background = element_rect(fill = "white", color = "black")) + ylim(c(0,1))
      
      ggsave(filename = paste0("yearseffect_",S_current,".jpeg"), y_plot)
    }
    
    # COVARIATE EFFECT ----------
    
    {
      years <- modelResults_MCMC$dataCharacteristics$Years
      Y <- length(years)
      
      X_tilde <- modelResults_MCMC$dataCharacteristics$X_tilde
      if(usingSpatial){
        X_centers <- nrow(X_tilde)  
      } else {
        X_centers <- 0
      }
      
      numTimeSpaceCov <- modelResults_MCMC$dataCharacteristics$numTimeSpaceCov
      
      namesCovariates <- modelResults_MCMC$dataCharacteristics$nameVariables_psi
      ncov_psi <- length(namesCovariates)
      
      if(usingSpatial){
        namesCovariates_xt <- c("X_T","Y_T")
      } else {
        namesCovariates_xt <- c()
      }
      
      namesCovariates_all <- c(namesCovariates_xt, namesCovariates)
      ncov_psi_all <- ncov_psi + numTimeSpaceCov
      
      # mcmc
      {
        beta_psi_output <- modelResults_MCMC$modelOutput$beta_psi_output
        
        beta_psi_output <- apply(beta_psi_output, 3, c)
        
        betaeffect_MCMC <- apply(beta_psi_output[,Y + X_centers + 1:ncov_psi_all,drop = F], 2, function(x){
          quantile(x, probs = c(0.025,0.5,0.975))
        })
      }
      
      # vb
      {
        mu_beta_psi <- modelResults_VB$modelOutput$mu_beta_psi
        Sigma_beta_psi <- modelResults_VB$modelOutput$Sigma_beta_psi
        
        betaeffect_VB <- sapply(1:(ncov_psi + numTimeSpaceCov), function(k){
          mu_beta_psi[Y + X_centers + k] + c(-1,1,0) * 1.96 * sqrt(Sigma_beta_psi[Y + X_centers + k,
                                                                                  Y + X_centers + k])
        })   
      }
      
      beta_effect_MCMC_all[,idx_s] <- betaeffect_MCMC
      beta_effect_VB_all[,idx_s] <- betaeffect_VB
      
      
    }
    
  }
  
}

# plots --------------

beta_psi_true_all <- beta_psi_true
if(usingSpatial){
  beta_psi_true_all <- c(beta_psi_tsp_true, beta_psi_true_all)
}

beta_MCMC_df <- as.data.frame(t(beta_effect_MCMC_all))
beta_MCMC_df$S <- factor(S_grid)

beta_VB_df <- as.data.frame(t(beta_effect_VB_all))
beta_VB_df$S <- factor(S_grid)

library(reshape2)

# beta_effect_MCMC_all_long_CI25 <- melt(beta_effect_MCMC_all[1,])
# beta_effect_MCMC_all_long_CI5 <- melt(beta_effect_MCMC_all[2,])
# beta_effect_MCMC_all_long_CI95 <- melt(beta_effect_MCMC_all[3,])

(occcov_plot <- ggplot() + 
  # geom_point(data = NULL, aes(x = namesCovariates_all,#colnames(betaeffect),
  #                             y = betaeffect[2,])) + 
  geom_errorbar(data = beta_MCMC_df, 
                aes(x = S,#colnames(betaeffect),
                                 ymin = V1,
                                 ymax = V3), alpha = .7) +
  geom_errorbar(data = beta_VB_df, aes(x = S,#colnames(betaeffect),
                                 ymin = V1,
                                 ymax = V3), size = 2) +
    scale_x_discrete() + 
  geom_hline(data = NULL, aes(yintercept = beta_psi_true), color = "red", size = 1) + 
  ylab("Effect") + scale_x_discrete(name = "Sites") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 16, face = "bold"),
        panel.grid.major = element_line(colour="grey", size=0.015),
        panel.background = element_rect(fill = "white", color = "black")))


ggsave(occcov_plot, file = "covariatesplot.jpeg")
