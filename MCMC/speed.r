library(microbenchmark)


microbenchmark({
   sample_z_cpp(psi, p, as.matrix(k_s[,3:4]))
},{
  XbetaY <- beta_psi[X_y_index]
  Xbeta <- XbetaY
  if(usingSpatial){
    if(spatialApprox == "SoD"){
      Xbetas <- beta_psi[Y + X_s_index]
    } else if(spatialApprox == "SoR"){
      Xbetas <- X_psi[,Y + seq_len(X_centers)] %*% beta_psi[Y + seq_len(X_centers)]
    }
    Xbeta <- Xbeta + Xbetas
  } else {
    Xbetas <- NULL
  }
  if((ncov_psi + numTimeSpaceCov) > 0) {
    Xbeta_cov <- X_psi[,-(1:(Y + X_centers)),drop = F] %*% beta_psi[-(1:(Y + X_centers))]
    Xbeta <- Xbeta + Xbeta_cov
  } else {
    Xbeta_cov <- rep(0, length(k))
  }
},{
  k <- z - .5
  n <- rep(1, length(k))
  # 
  # 
  #   list_beta_psi <- sampleBetaPsi_SoR(beta_psi, eps_s, X_psi, Xbeta, b_psi, inv_B_psi,
  #                                      n, k, Y, X_centers, ncov_psi, numTimeSpaceCov,
  #                                      X_y_index)
  # 
  Xbetaas <- Xbeta + eps_s
  Omega = samplePGvariables(Xbetaas, n)
},{
  XtOmegaX = matrixProductXtOmegaX_SoR(X_psi, Y, X_centers, numTimeSpaceCov + ncov_psi, 
                                   Omega, X_y_index, X_s_index)
},{
  knew = k - Omega * eps_s
  
  betaout = sample_beta_cpp_fast_sparse_SoR(X_psi, inv_B_psi, b_psi, knew, XtOmegaX,
                                                 X_y_index, X_s_index, Y, X_centers,
                                            numTimeSpaceCov + ncov_psi)
},{
  list_p <- update_p(beta_p, Occs, z_all, X_p, b_p, inv_B_p, Xbetap, 
                     ncov_p, p_intercepts, X_y_index_p, usingYearDetProb)
})




