load(here("Data/Simulated","simdata.rda"))

# years effect
{
  
  years <- modelResults_MCMC$dataCharacteristics$Years
  Y <- length(years)
  X_tilde <- modelResults_MCMC$dataCharacteristics$X_tilde
  if(usingSpatial){
    X_centers <- nrow(X_tilde)  
  } else {
    X_centers <- 0
  }
  
  beta_psi_output <- modelResults_MCMC$modelOutput$beta_psi_output
  beta_psi_output <- apply(beta_psi_output, 3, c)
  # beta_psi_output <- beta_psi_output[1,,]
  if(!usingSpatial){
    # a_s_unique_output <- modelResults$a_s_unique_output
    # a_s_unique_output <- apply(a_s_unique_output, 3, c)
  }
  
  CI_yearseffect <- sapply(1:Y, function(j) {
    
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
  
  prob_true <- logit(mu_psi_true + b_t_true #+ 
                     # beta_psi_t_true * sort(unique(modelResults$dataCharacteristics$X_psi_yearcov)) +
                     # mean(a_s_site_true)
  )
  
  ggplot(data = NULL, aes(x = years,
                          y = CI_yearseffect[3,],
                          ymin = CI_yearseffect[1,],
                          ymax = CI_yearseffect[2,])) + geom_point() + geom_errorbar() + geom_line() + 
    geom_point(data = NULL, aes(x = years, y = prob_true), size = 2, color = "red") +
    ylab("Occupancy probability") + scale_x_continuous(name = "Year", breaks = years) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black")) + ylim(c(0,1))
}

# spatial site patches
{
  X_tilde <- modelResults_MCMC$dataCharacteristics$X_tilde
  X_centers <- nrow(X_tilde)
  
  years <- modelResults_MCMC$dataCharacteristics$Years
  Y <- length(years)
  
  eps_unique_output <- modelResults_MCMC$modelOutput$eps_unique_output
  beta_psi_output <- modelResults_MCMC$modelOutput$beta_psi_output
  
  eps_unique_output <- apply(eps_unique_output, 2, c)
  beta_psi_output <- apply(beta_psi_output, 3, c)
  
  a_s_star <- sapply(1:nrow(X_tilde), function(i){
    mean(beta_psi_output[,Y + i] )#+ eps_unique_output[i])
    # )
  }) 
  
  ggplot() + geom_point(data = NULL, aes(x = X_tilde[,1],
                                         y = X_tilde[,2], color = a_s_star), 
                        size = 3, shape = 15) +
    labs(color="") + 
    scale_color_gradientn(colors = c("white","black")) +
    xlab("X") + ylab("Y") + theme(legend.title = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          # legend.position = "none",
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          panel.background = element_rect(fill = "white", color = "black"))
  
  # true
  
  ggplot() + 
    geom_point(data = NULL, aes(x = X[,1],
                                y = X[,2], alpha = a_s_site_true), size = 1) +
    xlab("X") + ylab("Y") + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          # legend.position = "none",
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          # panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
  
  if(usingSpatial){
    ggplot() +
      geom_point(data = NULL, aes(x = X[,1],
                                  y = X[,2], alpha = a_s_star), size = 1,
                 color = "black", shape = 15) +
      xlab("X") + ylab("Y") + theme(legend.title = element_blank()) +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            # legend.position = "none",
            axis.title = element_text(size = 20, face = "bold"),
            axis.text = element_text(size = 13, face = "bold", angle = 90),
            panel.background = element_rect(fill = "white", color = "black"))
  } else {
    estimatedPlot <- ggplot() + geom_point(data = NULL, aes(x = X_tilde_as[,2],
                                                            y = X_tilde_as[,3], alpha = a_s_mean), color = "black", size = 3) +
      xlab("X") + ylab("Y") + theme(legend.title = element_blank()) 
    
  }
  
  ggplot() + 
    geom_point(data = NULL, aes(x = X[,1],
                                y = X[,2], alpha = a_s_star - a_s_site_true), size = 1,
               color = "black") +
    xlab("X") + ylab("Y") + theme(legend.title = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          # legend.position = "none",
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          # panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
  
  mean(abs(a_s_site_true))
  mean(abs(a_s_star - a_s_site_true))
  
  qplot(a_s_site_true)
  qplot(a_s_star - a_s_site_true)
  
}

# spatial site effect + individual random effects
{
  X_tilde <- modelResults$dataCharacteristics$X_tilde
  X_centers <- nrow(X_tilde)
  
  years <- modelResults$dataCharacteristics$Years
  Y <- length(years)
  
  eps_unique_output <- modelResults$modelOutput$eps_unique_output
  beta_psi_output <- modelResults$modelOutput$beta_psi_output
  
  eps_unique_output <- apply(eps_unique_output, 2, c)
  beta_psi_output <- apply(beta_psi_output, 3, c)
  
  # X_sp <- data$X1
  # Y_sp <- data$X2
  # uniqueIndexesSite <- which(!duplicated(data$Site))
  # X_sp_unique <- X_sp[uniqueIndexesSite]
  # Y_sp_unique <- Y_sp[uniqueIndexesSite]
  # XY_sp_unique <- cbind(X_sp_unique, Y_sp_unique)
  # X_tilde_star <- XY_sp_unique
  
  #
  
  nSites <- nrow(X)
  idxSite <- findClosestPoint(X, X_tilde)
  a_s_star <- sapply(1:nSites, function(i){
    mean(beta_psi_output[,Y + idxSite[i]] + eps_unique_output[i])
    # )
  }) 
  
  ggplot() + geom_point(data = NULL, aes(x = X[,1],
                                         y = X[,2], alpha = a_s_site_true), size = 1,
                        color = "black") +
    xlab("X") + ylab("Y") + theme(legend.title = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          # legend.position = "none",
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          panel.background = element_rect(fill = "white", color = "black"))
  
  if(usingSpatial){
    ggplot() +
      geom_point(data = NULL, aes(x = X[,1],
                                  y = X[,2], alpha = a_s_star), size = 1,
                 color = "black", shape = 15) +
      xlab("X") + ylab("Y") + theme(legend.title = element_blank()) +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            # legend.position = "none",
            axis.title = element_text(size = 20, face = "bold"),
            axis.text = element_text(size = 13, face = "bold", angle = 90),
            panel.background = element_rect(fill = "white", color = "black"))
  } else {
    estimatedPlot <- ggplot() + geom_point(data = NULL, aes(x = X_tilde_as[,2],
                                                            y = X_tilde_as[,3], alpha = a_s_mean), color = "black", size = 3) +
      xlab("X") + ylab("Y") + theme(legend.title = element_blank()) 
    
  }
  
  ggplot() + 
    geom_point(data = NULL, aes(x = X[,1],
                                y = X[,2], alpha = a_s_star - a_s_site_true), size = 1,
               color = "black") +
    xlab("X") + ylab("Y") + theme(legend.title = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          # legend.position = "none",
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          # panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
  
  mean(abs(a_s_site_true))
  mean(abs(a_s_star - a_s_site_true))
  
  qplot(a_s_site_true)
  qplot(a_s_star - a_s_site_true)
  
}

# spatial site effect 
{
  gridStep <- modelResults_MCMC$dataCharacteristics$gridStep
  
  X_tilde <- modelResults_MCMC$dataCharacteristics$X_tilde
  X_tilde_as <- modelResults_MCMC$dataCharacteristics$X_tilde_as
  X_centers <- nrow(X_tilde)
  
  years <- modelResults_MCMC$dataCharacteristics$Years
  Y <- length(years)
  
  beta_psi_output <- modelResults_MCMC$modelOutput$beta_psi_output
  beta_psi_output <- apply(beta_psi_output, 3, c)
  
  X_sp <- data$X1
  Y_sp <- data$X2
  uniqueIndexesSite <- which(!duplicated(data$Site))
  X_sp_unique <- X_sp[uniqueIndexesSite]
  Y_sp_unique <- Y_sp[uniqueIndexesSite]
  XY_sp_unique <- cbind(X_sp_unique, Y_sp_unique)
  X_tilde_star <- XY_sp_unique
  
  #
  
  if(usingSpatial){
    mean_siteeffect <- sapply(1:X_centers, function(i){
      mean(beta_psi_output[,Y + i])
    })  
  } else {
    a_s_mean <- apply(modelResults$a_s_unique_output[1,,], 2, mean)
  }
  
  a_s_star <- mean_siteeffect[findClosestPoint(X, X_tilde)]
  
  #
  
  ggplot() + geom_point(data = NULL, aes(x = X[,1],
                                         y = X[,2], alpha = a_s_site_true), size = 1) +
    xlab("X") + ylab("Y") + theme(legend.title = element_blank()) +
    xlab("X") + ylab("Y") + #theme(legend.title = element_blank())
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.position = "none",
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          panel.background = element_rect(fill = "white", color = "black"))
  
  if(usingSpatial){
    ggplot() +
      geom_point(data = NULL, aes(x = X[,1],
                                  y = X[,2], alpha = a_s_star), size = 1) +
      xlab("X") + ylab("Y") + #theme(legend.title = element_blank())
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            legend.position = "none",
            axis.title = element_text(size = 20, face = "bold"),
            axis.text = element_text(size = 13, face = "bold", angle = 90),
            # panel.grid.major = element_line(colour="grey", size=0.015),
            panel.background = element_rect(fill = "white", color = "black"))
  } else {
    ggplot() + geom_point(data = NULL, aes(x = X_tilde_as[,2],
                                           y = X_tilde_as[,3], alpha = a_s_mean), color = "black", size = 3) +
      xlab("X") + ylab("Y") + theme(legend.title = element_blank()) 
    
  }
  
  ggplot() + 
    geom_point(data = NULL, aes(x = X_tilde_star[,1],
                                y = X_tilde_star[,2], alpha = a_s_star - a_s_site_true), 
               size = 1) +
    xlab("X") + ylab("Y") + #theme(legend.title = element_blank()) 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.position = "none",
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          # panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
  
  qplot(a_s_star - a_s_site_true)
  
  qplot(a_s_site_true)
  
}

# check individual effects
{
  gridStep <- modelResults$dataCharacteristics$gridStep
  
  X_tilde <- modelResults$dataCharacteristics$X_tilde
  X_tilde_as <- modelResults$dataCharacteristics$X_tilde_as
  X_centers <- nrow(X_tilde)
  
  years <- modelResults$dataCharacteristics$Years
  Y <- length(years)
  
  beta_psi_output <- modelResults$beta_psi_output
  
  beta_psi_output <- apply(beta_psi_output, 3, c)
  
  X_sp <- data$X1
  Y_sp <- data$X2
  uniqueIndexesSite <- which(!duplicated(data$Site))
  X_sp_unique <- X_sp[uniqueIndexesSite]
  Y_sp_unique <- Y_sp[uniqueIndexesSite]
  XY_sp_unique <- cbind(X_sp_unique, Y_sp_unique)
  X_tilde_star <- XY_sp_unique
  
  CI_siteeffect <- t(sapply(1:X_centers, function(i){
    quantile(beta_psi_output[,Y + i], probs = c(0.025,0.975))
  })  )
  
  idxPoints <- 501:600
  
  # X[i,]
  # X_tilde_star[i,]
  
  a_s_star <- CI_siteeffect[findClosestPoint(X[idxPoints,,drop = F], X_tilde),]
  
  ggplot(data = NULL, aes(x = idxPoints, y = a_s_site_true[idxPoints], 
                          ymin = a_s_star[,1],
                          ymax = a_s_star[,2])) + 
    geom_point(color = "red") + geom_errorbar() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          # legend.position = "none",
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          # panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
}

# covariates psi 
{
  X_tilde <- modelResults$dataCharacteristics$X_tilde
  if(usingSpatial){
    X_centers <- nrow(X_tilde)  
  } else {
    X_centers <- 0
  }
  
  years <- modelResults$dataCharacteristics$Years
  Y <- length(years)
  
  numTimeSpaceCov <- modelResults$dataCharacteristics$numTimeSpaceCov
  
  namesCovariates <- modelResults$dataCharacteristics$nameVariables_psi
  ncov_psi <- length(namesCovariates)
  
  if(usingSpatial){
    namesCovariates_xt <- c("X_T","Y_T")
  } else {
    namesCovariates_xt <- c()
  }
  
  namesCovariates_all <- c(namesCovariates_xt, namesCovariates)
  ncov_psi_all <- ncov_psi + numTimeSpaceCov
  
  beta_psi_output <- modelResults$modelOutput$beta_psi_output
  
  beta_psi_output <- apply(beta_psi_output, 3, c)
  
  betaeffect <- apply(beta_psi_output[,Y + X_centers + 1:ncov_psi_all,drop = F], 2, function(x){
    quantile(x, probs = c(0.025,0.5,0.975))
  })
  
  beta_psi_true_all <- beta_psi_true
  if(usingSpatial){
    beta_psi_true_all <- c(beta_psi_tsp_true, beta_psi_true_all)
  }
  
  betaeffect <- rbind(betaeffect, beta_psi_true_all)
  colnames(betaeffect) <- namesCovariates_all
  
  ggplot(data = NULL, aes(x = colnames(betaeffect),
                          y = betaeffect[2,],
                          ymin = betaeffect[1,],
                          ymax = betaeffect[3,])) + geom_point() + geom_errorbar() +
    geom_point(data = NULL, aes(x = colnames(betaeffect),
                                y = betaeffect[4,]), color = "red", size = 4) + 
    ylab("Effect") + scale_x_discrete(name = "Covariates") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 16, face = "bold"),
          panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
}

# capture probability
{
  beta_p_output <- modelResults$modelOutput$beta_p_output
  
  beta_p_output <- apply(beta_p_output, 3, c)
  
  peffect <- quantile(logit(beta_p_output[,1]), probs = c(0.025,0.5,0.975))
  
  ggplot(data = NULL, aes(x = "",
                          y = peffect[2],
                          ymin = peffect[1],
                          ymax = peffect[3])) + geom_point() + geom_errorbar() +
    geom_point(data = NULL, aes(x = "", y = logit(mu_p_true)), size = 5, color = "red") +
    ylab("Value") + scale_x_discrete(name = "Capture probability") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 16, face = "bold"),
          panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
}

# covariates p
{
  
  namesCovariates <- modelResults$dataCharacteristics$nameVariables_p
  ncov_p <- length(namesCovariates)
  
  beta_p_output <- modelResults$modelOutput$beta_p_output
  
  beta_p_output <- apply(beta_p_output, 3, c)
  
  betaeffect <- apply(beta_p_output[,-1,drop = F], 2, function(x){
    quantile(x, probs = c(0.025,0.5,0.975))
  })
  
  betaeffect <- rbind(betaeffect, beta_p_true)
  
  ggplot(data = NULL, aes(x = namesCovariates,
                          y = betaeffect[2,],
                          ymin = betaeffect[1,],
                          ymax = betaeffect[3,])) + geom_point() + geom_errorbar() +
    geom_point(data = NULL, aes(x = namesCovariates, y= betaeffect[4,]), color = "red", size = 4) +
    ylab("Effect") + scale_x_discrete(name = "Covariates") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 16, face = "bold"),
          panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
}

# -------


# spatial site effect
{
  gridStep <- modelResults$dataCharacteristics$gridStep
  
  X_tilde <- modelResults$dataCharacteristics$X_tilde
  X_tilde_as <- modelResults$dataCharacteristics$X_tilde_as
  X_centers <- nrow(X_tilde)
  
  years <- modelResults$dataCharacteristics$Years
  Y <- length(years)
  
  beta_psi_output <- modelResults$beta_psi_output
  
  beta_psi_output <- apply(beta_psi_output, 3, c)
  
  if(usingSpatial){
    # CI_siteeffect <- apply(beta_psi_output[,Y + 1:X_centers], 2, function(x){
    # mean(x)
    # })  
    CI_siteeffect <- sapply(1:X_centers, function(i){
      mean(beta_psi_output[,Y + i] #+ 
           # beta_psi_output[,Y + X_centers + 1] * X_tilde[i,1] +
           # beta_psi_output[,Y + X_centers + 2] * X_tilde[i,2]
      )
    })  
  } else {
    a_s_mean <- apply(modelResults$a_s_unique_output[1,,], 2, mean)
  }
  
  ggplot() + geom_point(data = NULL, aes(x = X[,1],
                                         y = X[,2], alpha = a_s_site_true), size = 3) +
    xlab("X") + ylab("Y") + theme(legend.title = element_blank()) 
  # ggplot() + geom_point(data = NULL, aes(x = X[,1],
  #                                        y = X[,2], color = a_s_site_true)) +
  #   xlab("X") + ylab("Y") + theme(legend.title = element_blank()) 
  
  if(usingSpatial){
    ggplot() + geom_point(data = NULL, aes(x = X_tilde[,1],
                                           y = X_tilde[,2], alpha = CI_siteeffect), 
                          color = "black", size = 3.5 / gridStep, shape = 15)   +
      # geom_point(data = NULL, aes(x = X[,1],
      # y = X[,2], alpha = a_s_site_true), color = "red", size = 3) +
      xlab("X") + ylab("Y") + #theme(legend.title = element_blank()) 
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.title = element_text(size = 20, face = "bold"),
            axis.text = element_text(size = 13, face = "bold", angle = 90),
            # panel.grid.major = element_line(colour="grey", size=0.015),
            panel.background = element_rect(fill = "white", color = "black"))
  } else {
    ggplot() + geom_point(data = NULL, aes(x = X_tilde_as[,2],
                                           y = X_tilde_as[,3], alpha = a_s_mean), color = "black", size = 3) +
      xlab("X") + ylab("Y") + theme(legend.title = element_blank()) 
    
  }
  
}