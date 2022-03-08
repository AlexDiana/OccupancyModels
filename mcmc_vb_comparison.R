load(here("Data/Simulated","simdata.rda"))

setwd(here())

# YEAR EFFECT ---------

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
  
  ggsave(filename = paste0("yearseffect_",1000,".jpeg"), y_plot)
  
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
    
    # betaeffect <- rbind(betaeffect, beta_psi_true_all)
    # colnames(betaeffect) <- namesCovariates_all
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
  
  beta_psi_true_all <- beta_psi_true
  if(usingSpatial){
    beta_psi_true_all <- c(beta_psi_tsp_true, beta_psi_true_all)
  }
  
  
  occcov_plot <- ggplot() + 
    # geom_point(data = NULL, aes(x = namesCovariates_all,#colnames(betaeffect),
    #                             y = betaeffect[2,])) + 
    geom_errorbar(data = NULL, aes(x = namesCovariates_all,#colnames(betaeffect),
                                   ymin = betaeffect_MCMC[1,],
                                   ymax = betaeffect_MCMC[3,]), alpha = .7) +
    geom_errorbar(data = NULL, aes(x = namesCovariates_all,#colnames(betaeffect),
                                   ymin = betaeffect_VB[1,],
                                   ymax = betaeffect_VB[3,]), size = 2) +
    geom_point(data = NULL, aes(x = namesCovariates_all,
                                y = beta_psi_true), color = "red", size = 4) + 
    ylab("Effect") + scale_x_discrete(name = "Covariates") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 16, face = "bold"),
          panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
  
  ggsave(filename = paste0("occ_covariates_",1000,".jpeg"), occcov_plot)
  
}

# DETECTION PROBABILITY ---------

# MCMC
{
  beta_p_output <- modelResults_MCMC$modelOutput$beta_p_output
  
  beta_p_output <- apply(beta_p_output, 3, c)
  
  if(usingYearDetProb){
    peffect_MCMC <- apply(beta_p_output[,1:Y], 2, function(x) {
      quantile(logit(x), probs = c(0.025,0.5,0.975))
    })
  } else {
    peffect_MCMC <- quantile(logit(beta_p_output[,1]), probs = c(0.025,0.5,0.975))
  }
  
  
}

# VB 
{
  mu_beta_p <- modelResults_VB$modelOutput$mu_beta_p
  Sigma_beta_p <- modelResults_VB$modelOutput$Sigma_beta_p
  
  if(usingYearDetProb){
    peffect_VB <- sapply(1:Y, function(j){
      logit(mu_beta_p[j] + c(-1,1,0) * 1.96 * sqrt(Sigma_beta_p[j,j]))
    })
  } else {
    peffect_VB <- logit(mu_beta_p[1] + c(-1,1,0) * 1.96 * sqrt(Sigma_beta_p[1,1]))
  }
}

detprob_plot <- ggplot() + 
  geom_errorbar(data = NULL, aes(x = "",
                                 y = peffect_MCMC[2],
                                 ymin = peffect_MCMC[1],
                                 ymax = peffect_MCMC[3]), alpha = .7) +
  geom_errorbar(data = NULL, aes(x = "",
                                 y = peffect_VB[2],
                                 ymin = peffect_VB[1],
                                 ymax = peffect_VB[3]), size = 2) +
  geom_point(data = NULL, aes(x = "", y = logit(mu_p_true)), size = 5, color = "red") +
  ylab("Value") + scale_x_discrete(name = "Capture probability") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 16, face = "bold"),
        panel.grid.major = element_line(colour="grey", size=0.015),
        panel.background = element_rect(fill = "white", color = "black"))

ggsave(filename = paste0("det_prob_",1000,".jpeg"), detprob_plot)

