

tanh <- function(x) (exp(2*x) - 1) / (exp(2*x) + 1)

update_p <- function(mu_beta_p, Sigma_beta_p, 
                     y, gamma_z_all, X_p, b_p, inv_B_p, Xbeta_p, 
                     ncov_p, p_intercepts, X_y_index_p, usingYearDetProb){
  
  Sigma_beta_p_0 <- Sigma_beta_p  
  mu_beta_p_0 <- mu_beta_p
  # update pg variational params
  {
    phi_p <- sapply(1:length(y), function(i){
      -.5 * (t(X_p[i,]) %*% Sigma_beta_p_0 %*% X_p[i,] + (t(X_p[i,]) %*% mu_beta_p_0)^2)
    })
    
    csi_p <- sqrt(-2 * phi_p)
  }
  
  # update beta
  {
    diagZeta <- gamma_z_all * .5 *tanh(.5 * csi_p) / csi_p
    X_pZetaX_p <- diagMatrixProd(t(X_p), diagZeta) %*% X_p

    Lambda_2 <- - .5 * (inv_B_p + X_pZetaX_p) 

    Sigma_beta_p <- solve(- 2 * Lambda_2)
    
    Lambda_1 <- t(X_p) %*% (gamma_z_all * (y - .5)) + inv_B_p %*% b_p
    mu_beta_p <- Sigma_beta_p %*% Lambda_1
  }
  
  #
  list_new <- update_p(mu_beta_p_0, Sigma_beta_p_0, y, gamma_z_all, X_p, b_p, 
                         inv_B_p, mu_beta_p, ncov_p, p_intercepts, X_y_index_p, usingYearDetProb)
  
  list_new$mu_beta_p
  
  # k <- y - .5
  # n <- rep(1, length(k))
  # 
  # k <- k[z_all==1]
  # n <- n[z_all==1]
  # X_p_present <- X_p[z_all == 1,,drop = FALSE]
  # X_y_index_p_present <- X_y_index_p[z_all == 1]
  # 
  # Xbeta <- Xbeta_p[z_all == 1]
  # 
  # # list_beta_p <- sample_beta_omega_cpp_parallel(beta_p, X_p_present, b_p, B_p, n, k)
  # # list_beta_p <- sample_beta_omega_cpp(beta_p, X_p_present, b_p, inv_B_p, n, k)
  # # beta_p <- list_beta_p$beta
  # beta_p <- sampleBetaP(beta_p, X_p_present, Xbeta, b_p, inv_B_p,
  #                         ncov_p, p_intercepts, X_y_index_p_present, n, k)
  # 
  # Xbeta <- X_p[,-(1:p_intercepts),drop = F] %*% beta_p[-(1:p_intercepts)]
  # XbetaY <- beta_p[X_y_index_p]
  # Xbeta <- Xbeta + XbetaY
  # 
  # p <- as.vector(logit(Xbeta))
  
  # p <- as.vector(logit(X_p %*% beta_p))
  
  list("mu_beta_p" = mu_beta_p, "Sigma_beta_p" = Sigma_beta_p)
}
