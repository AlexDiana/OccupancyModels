
update_psi <- function(beta_psi, X_psi_cov, Xbeta, 
                       b_psi, inv_B_psi, 
                       z, k_s, sites, Y, X_centers, ncov_psi, 
                       X_y_index, X_sc_index, X_s_sor,
                       numTimeSpaceCov,
                       usingSpatial, spatialApprox, eps_s, sigma_eps){
  
  k <- z - .5
  n <- rep(1, length(k))
  
  # list_beta_psi <- sampler_beta(beta_psi, a_s, X_psi, b_psi, inv_B_psi,
  #                                  n, k, Y, X_centers, ncov_psi, numTimeSpaceCov,
  #                                  X_y_index, X_s_index)
  # list_beta_psi <- sampler_beta_parallel(beta_psi, eps_s, X_psi, b_psi, inv_B_psi,
                                         # n, k, Y, X_centers, ncov_psi, numTimeSpaceCov,
                                         # X_y_index, X_s_index)
  if(!usingSpatial | (usingSpatial & spatialApprox == "SoD")){
    list_beta_psi <- sampleBetaPsi(beta_psi, eps_s, X_psi_cov, Xbeta, b_psi, inv_B_psi,
                                   n, k, Y, X_centers, ncov_psi, numTimeSpaceCov,
                                   X_y_index,  X_sc_index)
  } else {
    list_beta_psi <- sampleBetaPsi_SoR(beta_psi, eps_s, X_psi_cov, Xbeta, b_psi, inv_B_psi,
                                       n, k, Y, X_centers, ncov_psi, numTimeSpaceCov,
                                       X_y_index, X_sc_index, X_s_sor)
  }
    
  beta_psi <- list_beta_psi$beta
  Omega <- list_beta_psi$Omega
  
  # compute X * beta and store the results for later use
  
  XbetaY <- beta_psi[X_y_index]
  Xbeta <- XbetaY
  if(usingSpatial){
    if(spatialApprox == "SoD"){
      Xbetas <- beta_psi[Y + X_sc_index]
    } else if(spatialApprox == "SoR"){
      # Xbetas <- X_psi[,Y + seq_len(X_centers)] %*% beta_psi[Y + seq_len(X_centers)]
      Xbetas <- X_s_sor %*% beta_psi[Y + seq_len(X_centers)]
    }
    Xbeta <- Xbeta + Xbetas
  } else {
    Xbetas <- NULL
  }
  if((ncov_psi + numTimeSpaceCov) > 0) {
    Xbeta_cov <- X_psi_cov %*% beta_psi[-(1:(Y + X_centers))]
    Xbeta <- Xbeta + Xbeta_cov
  } else {
    Xbeta_cov <- rep(0, length(k))
  }
  
  # update random effects
  
  # eps_s <- sample_eps_cpp(k_s$Site, sort(unique(k_s$Site)),
  #                         Xbeta, k,
  #                         z, Omega, sigma_eps)

  psi <- as.vector(logit(Xbeta + eps_s))
  
  list("beta" = beta_psi, "psi" = psi, "Omega" = Omega, "eps_s" = eps_s,
       "XbetaY" = XbetaY, "Xbetas" = Xbetas, "Xbeta_cov" = Xbeta_cov, "Xbeta" = Xbeta)
}
