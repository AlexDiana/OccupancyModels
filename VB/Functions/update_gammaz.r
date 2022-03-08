update_gamma_z <- function(gamma_z, k_s, mu_eps, 
                           mu_beta_psi, X_psi,
                           X_p, mu_beta_p){
  
  l <- 0
  
  for (j in 1:length(gamma_z)) {
    
    if(k_s[j,3] == 0){
      
      firsTerm <- t(X_psi[j,]) %*% mu_beta_psi + mu_eps[k_s[j,2]]
      secondTerm <- log(1 + exp(t(X_psi[j,]) %*% mu_beta_psi + mu_eps[k_s[j,2]]))
      # secondTerm <- mean(sapply(1:10, function(i){
      #   beta_psi <- mvrnorm(1, mu_beta_psi, Sigma_beta_psi)
      #   eps_j <- rnorm(1, mu_eps[k_s[j,2]], sd_eps[k_s[j,2]])
      #   log(1 + exp(t(X_psi[j,]) %*% beta_psi + eps_j))
      # }))
      
      meanlogphi <- firsTerm - secondTerm

      meanlog1mphi <- - log(1 + exp(t(X_psi[j,]) %*% mu_beta_psi + mu_eps[k_s[j,2]]))
      # meanlog1mphi <- mean(sapply(1:10, function(j){
      #   beta_psi <- mvrnorm(1, mu_beta_psi, Sigma_beta_psi)
      #   - log(1 + exp(t(X_psi[j,]) %*% mu_beta_psi + mu_eps[k_s[j,2]]))
      # })
                  
      which_kij <- l + 1:k_s[j,4]
      
      c_i <- sapply(which_kij, function(i){
        # mean(sapply(1:10,function(k){
        #   beta_p <- mvrnorm(1, mu_beta_p, Sigma_beta_p)
        #   - log(1 + exp(t(X_p[i,]) %*% beta_p))
        # }))
        - log(1 + exp(t(X_p[i,]) %*% mu_beta_p))
      })
      sum(c_i)
      
      a_star <- sum(c_i) + meanlogphi
      b_star <- meanlog1mphi
      c_star <- a_star - b_star
      
      gamma_z[j] <- exp(c_star) / (1 + exp(c_star))
      
    } else {
      
      gamma_z[j] <- 1
      
    }
    
    l <- l + k_s[j,4]
    
  }
  
  gamma_z
}

