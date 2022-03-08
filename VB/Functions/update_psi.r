

tanh <- function(x) (exp(2*x) - 1) / (exp(2*x) + 1)

mean_pg <- function(c){
  tanh(c / 2) / (2 * c)
}

update_psi <- function(mu_beta_psi, Sigma_beta_psi,
                       X_psi, Xbeta, 
                       b_psi, inv_B_psi, 
                       gamma_z, k_s, sites, Y, X_centers, ncov_psi, 
                       X_y_index, X_s_index, numTimeSpaceCov,
                       usingSpatial, 
                       mu_eps, sd_eps, sigma_eps,
                       a_l_T, b_l_T, a_sigma_T, b_sigma_T,
                       a_l_s, b_l_s, a_sigma_s, b_sigma_s){
  
  
  # update pg variational params
  {
    phi_psi <- sapply(1:nrow(k_s), function(j){
      # -.5 * (t(X_psi[i,]) %*% Sigma_beta_psi %*% X_psi[i,] + (t(X_psi[i,]) %*% mu_beta_psi)^2) 
      Xit_mu <- t(X_psi[j,]) %*% mu_beta_psi
      XiiSigmaXi <- t(XiSigmapsi(X_psi[j,], X_y_index[j], X_s_index[j], Sigma_beta_psi, 
                               Y, X_centers, numTimeSpaceCov + ncov_psi)) %*% X_psi[j,]
      # -.5 * (t(X_psi[j,]) %*% Sigma_beta_psi %*% X_psi[j,] + Xit_mu^2 +
               # 2 * Xit_mu * mu_eps[j] + sd_eps[j]^2 + mu_eps[j]^2)
      -.5 * (XiiSigmaXi + Xit_mu^2 +
               2 * Xit_mu * mu_eps[j] + sd_eps[j]^2 + mu_eps[j]^2)
    })
    
    csi_psi <- sqrt(-2 * phi_psi)
  }
  
  # update beta
  {
    # Zeta <- .5 * diag(tanh(.5 * csi_psi) / csi_psi)
    diagZeta <- .5 * tanh(.5 * csi_psi) / csi_psi
    X_psiZetaX_psi <- diagMatrixProd(t(X_psi), diagZeta) %*% X_psi
    
    X_psiZeta <- diagMatrixProd(t(X_psi), diagZeta)
    
    X_psiZeta <- t(X_psi[j,]) %*% Sigma_beta_psi
    
    X_psiZeta2 <- XiSigmapsi(X_psi[j,], X_y_index[j], X_s_index[j], Sigma_beta_psi, 
                             Y, X_centers, numTimeSpaceCov + ncov_psi)
    
    
    
    max(abs(X_psiZeta[1,] - X_psiZeta2[,1]))
    
    diagMatrixProd(t(X_psi), diagZeta)
    
    
    # sims_invB <- 5000
    # inv_B_psi_list <- lapply(1:sims_invB, function(i){
    #   l_T <- rgamma(1, a_l_T, b_l_T)
    #   sigma_T <- rgamma(1, a_sigma_T, b_sigma_T)
    #   l_S <- rgamma(1, a_l_s, b_l_s)
    #   sigma_S <- rgamma(1, a_sigma_s, b_sigma_s)
    #   
    #   inv_B_psi_new <- inv_B_psi
    #   
    #   K_l <- K(1:Y, 1:Y, sigma_T^2, l_T) + sigma_psi^2 + diag(exp(-8), nrow = Y)
    #   
    #   inverse_Kl <- FastGP::rcppeigen_invert_matrix(K_l)
    #   inverse_Kl[lower.tri(K_l)] <- t(inverse_Kl)[lower.tri(K_l)]
    #   
    #   inv_B_psi_new[1:Y, 1:Y] <- inverse_Kl
    #   
    #   if(usingSpatial){
    #     K_S <- K(X_centers, X_centers, sigma_S^2, l_S) + diag(exp(-8), nrow = Y)
    #     
    #     inverse_Kl <- FastGP::rcppeigen_invert_matrix(K_S)
    #     inverse_Kl[lower.tri(K_S)] <- t(inverse_Kl)[lower.tri(K_S)]
    #     
    #     inv_B_psi_new[Y + 1:X_centers, Y + 1:X_centers] <- inverse_Kl
    #   }
    #   
    #   inv_B_psi_new
    # })  
    # 
    # invB_psi_mean <- Reduce("+", inv_B_psi_list) / sims_invB
    
    # Lambda_2 <- - .5 * (inv_B_psi + X_psiZetaX_psi)
    # Lambda_2 <- - .5 * (invB_psi_mean + X_psiZetaX_psi)
    # Sigma_beta_psi <- solve(- 2 * Lambda_2)
    
    Sigma_beta_psi <- solve(inv_B_psi + X_psiZetaX_psi)
    
    y_tilde <- (gamma_z - .5) - mean_pg(csi_psi) * mu_eps
    Lambda_1 <- t(X_psi) %*% y_tilde + inv_B_psi %*% b_psi
    # Lambda_1 <- t(X_psi) %*% y_tilde + invB_psi_mean %*% b_psi
    (mu_beta_psi <- Sigma_beta_psi %*% Lambda_1)
  }
  
  # update eps
  {
    X_psi_mu_beta_psi <- X_psi %*% mu_beta_psi
    mu_sd_eps <- sapply(1:nrow(k_s), function(j){
      
      mean_xib <- X_psi_mu_beta_psi[j]
      b <- (gamma_z[j] - .5) - mean_pg(csi_psi[j]) * mean_xib
      a <- mean_pg(csi_psi[j]) + 1 / sigma_eps^2
      
      c( b / a, sqrt(1 / a))
    })
    
    mu_eps <- mu_sd_eps[1,]
    sd_eps <- mu_sd_eps[2,]
  }
  
  list("mu_beta_psi" = mu_beta_psi, "Sigma_beta_psi" = Sigma_beta_psi,
       "mu_eps" = mu_eps, "sd_eps" = sd_eps)
}
