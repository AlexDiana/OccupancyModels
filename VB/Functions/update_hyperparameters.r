logp_T_r <- function(y, l, sig, x_gp, sig2){
  K_l <- K(x_gp, x_gp, sig, l) + sig2^2 + .001 * diag(1, nrow = length(x_gp))
  rcpp_log_dmvnorm(K_l, rep(0, length(x_gp)), y, F) 
}

logp_s_r <- function(y, l, sig, x_gp, sig2){
  K_l <- K2(x_gp, x_gp, sig, l) + .001 * diag(1, nrow = nrow(x_gp))
  rcpp_log_dmvnorm(K_l, rep(0, nrow(x_gp)), y, F)
}

logp <- function(y, l, sig, x_gp, sig2){
  K_l <- K2(x_gp, x_gp, sig, l) + sig2^2 + .001 * diag(1, nrow = nrow(x_gp))
  rcpp_log_dmvnorm(K_l, rep(0, length(x_gp)), y, F) 
}

logp0 <- function(x, a_x, b_x){
  dgamma(x, a_x, b_x, log = T)
}

logq <- function(x, a, b){
  dgamma(x, a, b, log = T)
}

Grad_logq <- function(x, a, b){
  c(log(b) - digamma(a) + log(x),
    a / b - x)
}

G_lambda <- function(a, b){
  matrix(c(trigamma(a), - 1 / b,
           - 1 / b, a / b^2), 2, 2)
}

# NEW METHOD -----

logp0 <- function(x, a_x, b_x){
  dgamma(x, a_x, b_x, log = T)
}

logq <- function(x, a, b){
  dgamma(x, a, b, log = T)
}

Grad_logq <- function(x, a, b){
  c(log(b) - digamma(a) + log(x),
    a / b - x)
}

eval_ELBO_T <- function(l_T,
                      a_l_0, b_l_0,
                      a_sigma_T, b_sigma_T, 
                      y, x, sig, sims_grad){
  
  grad_sims <- sapply(1:sims_grad, function(s){
    # print(s)
    sig_T <- rgamma(1, a_sigma_T, b_sigma_T)
    
    logp_term <- logp_T(y,
                        l_T, 
                        sig_T^2, x, sig) 
    
    logp_term  + 
      dgamma(l_T, a_l_0, b_l_0, log = T)
  })
  
  mean(grad_sims)
  
}

eval_ELBO_S <- function(l_S,
                      a_l_0, b_l_0,
                      a_sigma_T, b_sigma_T, 
                      y, x, sig, sims_grad){
  
  grad_sims <- sapply(1:sims_grad, function(s){
    # print(s)
    sig_T <- rgamma(1, a_sigma_T, b_sigma_T)
    
    logp_term <- logp_s(y,
                        l_T, 
                        sig_T^2, x, sig) 
    
    logp_term  + 
      dgamma(l_S, a_l_0, b_l_0, log = T)
  })
  
  mean(grad_sims)
  
}


update_l_old <- function(y, a_l, b_l, 
                     a_sig, b_sig,
                     sig2, x_gp,
                     a_l_0, b_l_0,
                     sims_grad, rho){
  
  grad_sims <- sapply(1:sims_grad, function(s){
    # print(s)
    l_s <- rgamma(1, a_l_S, b_l_S)
    sig_s <- rgamma(1, a_sigma_S, b_sigma_S)
    h_i <- Grad_logq(l_s, a_l_S, b_l_S)
    d_i <- logp(y, l_s, sig_s, x_gp, sig2) +
      logp0(l_s, a_l_0, b_l_0) - logq(l_s, a_l, b_l)
    f_i <- h_i %*% d_i
    cbind(h_i, rep(d_i, length(f_i)), f_i)
  })
  
  a_i <- sum(sapply(1:2, function(d){
    cov(grad_sims[d + 4,], grad_sims[d,])
  })) /  sum(sapply(1:2, function(d){
    var(grad_sims[d,])
  }))
  
  grad_sims_new <- sapply(1:sims_grad, function(s){
    grad_sims[1:2, s] * (grad_sims[3,s] - a_i)   
  })
  
  grad_est <- apply(grad_sims_new, 1, mean)
 
  Glambda <- G_lambda(a_l, b_l)
  
  nat_grad <- solve(Glambda) %*% grad_est
  
  a_l <-  a_l + rho * nat_grad[1]
  b_l <-  b_l + rho * nat_grad[2]
  
  list("a_l" = a_l,
       "b_l" = b_l)
}


update_l_T <- function(y, l_T, 
                     a_sig, b_sig,
                     sig2, x_gp,
                     a_l_0, b_l_0,
                     sims_grad){
  
  pointsToSearch <- c(seq(l_T / 2, l_T, length.out = 5),
                      seq(l_T * 1.1, l_T * 10, length.out = 5))
  
  ( elbo_vals <- sapply(1:length(pointsToSearch), function(i){
    # print(i)
    
    eval_ELBO_T(pointsToSearch[i],
              a_l_0, b_l_0,
              a_sig, b_sig,
              y,  x_gp, sig2, sims_grad)
  }) )
  
  l_T <- pointsToSearch[which.max(elbo_vals)]
  
  l_T
}

update_l_sigma_T <- function(y, l_T, 
                     sigma_T,
                     sig2, x_gp,
                     a_l_T_0, b_l_T_0,
                     a_sigma_T_0, b_sigma_T_0){
  
  pointsToSearch_sigma <- c(seq(sigma_T / 2, sigma_T, length.out = 5),
                      seq(sigma_T, sigma_T * 10, length.out = 5 + 1)[-1])
  
  pointsToSearch_l <- c(seq(l_T / 5, l_T, length.out = 5),
                        seq(l_T, l_T * 10, length.out = 5 + 1)[-1])
  
  pointsToSearch <- expand.grid(pointsToSearch_sigma,
                                pointsToSearch_l)
  
  ( elbo_vals <- sapply(1:nrow(pointsToSearch), function(i){
    # print(i)
    
    sigma_T_current <- pointsToSearch[i,1]
    l_T_current <- pointsToSearch[i,2]
    
    logp_term <- logp_T(y,
                        l_T_current, 
                        sigma_T_current^2, x_gp, sig2) 
    
    logp_term  + 
      dgamma(l_T_current, a_l_T_0, b_l_T_0, log = T) +
      dgamma(sigma_T_current, a_sigma_T_0, b_sigma_T_0, log = T)
  
  }) )
  
  sigma_T <- pointsToSearch[which.max(elbo_vals), 1]
  l_T <- pointsToSearch[which.max(elbo_vals), 2]
  
  list("sigma_T" = sigma_T,
       "l_T" = l_T)
}

update_l_T_new <- function(y, l_T, 
                     a_sigma_T, b_sigma_T,
                     sig2, x_gp,
                     a_l_T_0, b_l_T_0,
                     a_sigma_T_0, b_sigma_T_0){
  
  sigmasq_T <- b_sigma_T / (a_sigma_T - 1)
  
  pointsToSearch <- c(seq(l_T / 5, l_T, length.out = 5),
                        seq(l_T, l_T * 10, length.out = 5 + 1)[-1])
  
  ( elbo_vals <- sapply(1:length(pointsToSearch), function(i){
    # print(i)
    
    l_T_current <- pointsToSearch[i]
    
    logp_term <- logp_T(y,
                        l_T_current, 
                        sigmasq_T, x_gp, sig2) 
    
    logp_term  + 
      dgamma(l_T_current, a_l_T_0, b_l_T_0, log = T) 
  
  }) )
  
  l_T <- pointsToSearch[which.max(elbo_vals)]
  
  l_T
}

update_l_S_new <- function(y, l_S, 
                     a_sigma_S, b_sigma_S,
                     sig2, x_gp,
                     a_l_S_0, b_l_S_0,
                     a_sigma_S_0, b_sigma_S_0){
  
  sigmasq_S <- b_sigma_S / (a_sigma_S - 1)
  
  pointsToSearch <- c(seq(l_S / 5, l_S, length.out = 3),
                        seq(l_S, l_S * 10, length.out = 3 + 1)[-1])
  
  ( elbo_vals <- sapply(1:length(pointsToSearch), function(i){
    # print(i)
    
    l_S_current <- pointsToSearch[i]
    
    logp_term <- logp_s(y,
                        l_S_current, 
                        sigmasq_S, x_gp, sig2) 
    
    logp_term  + 
      dgamma(l_S_current, a_l_S_0, b_l_S_0, log = T) 
  
  }) )
  
  l_S <- pointsToSearch[which.max(elbo_vals)]
  
  l_S
}

update_l_S_grid <- function(y, 
                            l_s_grid, diag_K_s_grid,
                            inv_K_s_grid, 
                            a_sigma_S, b_sigma_S,
                            a_l_S_0, b_l_S_0){
  
  sigmasq_S <- b_sigma_S / (a_sigma_S - 1)
  
  pointsToSearch <- l_s_grid
  
  ( elbo_vals <- sapply(1:length(pointsToSearch), function(i){
    # print(i)
    
    l_S_current <- pointsToSearch[i]
    
    # logp_term <- logp_s(y,
    #                     l_S_current, 
    #                     sigmasq_S, x_gp, sig2) 
    
    logp_term <- ((-length(y)/2) * log(2 * pi) -
                    sum(log(diag_K_s_grid[,i]) + log(sqrt(sigmasq_S))) -
                    (1/2) * (1 / sigmasq_S) * y %*% inv_K_s_grid[,,i] %*% y)
    
    logp_term  + 
      dgamma(l_S_current, a_l_S_0, b_l_S_0, log = T) 
  
  }) )
  
  l_S <- pointsToSearch[which.max(elbo_vals)]
  inv_Kl <- inv_K_s_grid[,,which.max(elbo_vals)]
    
  list("l_S" = l_S,
       "inv_Kl" = inv_Kl)
}

update_l_S <- function(y, l_S, 
                     a_sig, b_sig,
                     sig2, x_gp,
                     a_l_0, b_l_0,
                     sims_grad){
  
  pointsToSearch <- c(seq(l_S / 2, l_S, length.out = 5),
                      seq(l_S * 1.1, l_S * 10, length.out = 5))
  
  ( elbo_vals <- sapply(1:length(pointsToSearch), function(i){
    # print(i)
    
    eval_ELBO_S(pointsToSearch[i],
              a_l_0, b_l_0,
              a_sig, b_sig,
              y,  x_gp, sig2, sims_grad)
  }) )
  
  l_S <- pointsToSearch[which.max(elbo_vals)]
  
  l_S
}

update_sig <- function(y, a_l, b_l, 
                       a_sig, b_sig,
                       sig2, x_gp,
                       a_sig_0, b_sig_0,
                       sims_grad, rho){
  
  grad_sims <- sapply(1:sims_grad, function(s){
    # print(s)
    l_s <- rgamma(1, a_l, b_l)
    sig_s <- rgamma(1, a_sig, b_sig)
    h_i <- Grad_logq(sig_s, a_sig, b_sig)
    d_i <- logp(y, l_s, sig_s, x_gp, sig2) +
      logp0(sig_s, a_sig_0, b_sig_0) - logq(sig_s, a_sig, b_sig)
    f_i <- h_i %*% d_i
    cbind(h_i, rep(d_i, length(f_i)), f_i)
  })
  
  a_i <- sum(sapply(1:2, function(d){
    cov(grad_sims[d + 4,], grad_sims[d,])
  })) /  sum(sapply(1:2, function(d){
    var(grad_sims[d,])
  }))
  
  grad_sims_new <- sapply(1:sims_grad, function(s){
    grad_sims[1:2, s] * (grad_sims[3,s] - a_i)   
  })
  
  grad_est <- apply(grad_sims_new, 1, mean)
 
  Glambda <- G_lambda(a_sig, b_sig)
  
  nat_grad <- solve(Glambda) %*% grad_est
  
  a_sig <-  a_sig + rho * nat_grad[1]
  b_sig <-  b_sig + rho * nat_grad[2]
  
  list("a_sig" = a_sig,
       "b_sig" = b_sig)
}


update_l_sigma_S_old <- function(y, l_S, 
                             sigma_S,
                             sig2, x_gp,
                             a_l_S_0, b_l_S_0,
                             a_sigma_S_0, b_sigma_S_0){
  
  pointsToSearch_sigma <- c(seq(sigma_S / 2, sigma_S, length.out = 5),
                            seq(sigma_S * 1.1, sigma_S * 10, length.out = 5))
  
  pointsToSearch_l <- c(seq(l_S / 2, l_S, length.out = 5),
                        seq(l_S * 1.1, l_S * 10, length.out = 5))
  
  pointsToSearch <- expand.grid(pointsToSearch_sigma,
                                pointsToSearch_l)
  
  ( elbo_vals <- sapply(1:nrow(pointsToSearch), function(i){
    # print(i)
    
    sigma_S_current <- pointsToSearch[i,1]
    l_S_current <- pointsToSearch[i,2]
    
    logp_term <- logp_s_r(y,
                        l_S_current, 
                        sigma_S_current^2, x_gp, sig2) 
    
    logp_term  + 
      dgamma(l_S_current, a_l_S_0, b_l_S_0, log = T) +
      dgamma(sigma_S_current, a_sigma_S_0, b_sigma_S_0, log = T)
    
  }) )
  
  sigma_S <- pointsToSearch[which.max(elbo_vals), 1]
  l_S <- pointsToSearch[which.max(elbo_vals), 2]
  
  list("sigma_S" = sigma_S,
       "l_S" = l_S)
}

update_l_sigma_S <- function(y, l_S, 
                             sigma_S,
                             sig2, x_gp,
                             a_l_S_0, b_l_S_0,
                             a_sigma_S_0, b_sigma_S_0){
  
  pointsToSearch_sigma <- c(seq(sigma_S / 2, sigma_S, length.out = 3),
                            seq(sigma_S * 1.1, sigma_S * 10, length.out = 4))
  
  pointsToSearch_l <- c(seq(l_S / 2, l_S, length.out = 3),
                        seq(l_S * 1.1, l_S * 10, length.out = 4))
  
  pointsToSearch <- expand.grid(pointsToSearch_sigma,
                                pointsToSearch_l)
  
  elbo_vals <- computeElboVals(mu_beta_psi[Y + 1:X_centers], 
                                as.matrix(pointsToSearch), 0, X_tilde,
                                a_l_s_0, b_l_s_0, 
                                a_sigma_s_0, b_sigma_s_0)
  
  sigma_S <- pointsToSearch[which.max(elbo_vals), 1]
  l_S <- pointsToSearch[which.max(elbo_vals), 2]
  
  list("sigma_S" = sigma_S,
       "l_S" = l_S)
}


# OLD ONE -------

# update_l_T <- function(y, a_l, b_l, 
#                      a_sig, b_sig,
#                      sig2, x_gp,
#                      a_l_0, b_l_0,
#                      G_t,
#                      sims_grad, rho){
#   
#   grad_est <- apply(sapply(1:sims_grad, function(s){
#     l_s <- rgamma(1, a_l, b_l)
#     sig_s <- rgamma(1, a_sig, b_sig)
#     Grad_logq(l_s, a_l, b_l) %*% (logp_T(y, l_s, sig_s, x_gp, sig2) +
#                                     logp0(l_s, a_l_0, b_l_0) - logq(l_s, a_l, b_l))
#   }), 1, mean)
#   
#   # if(all(diag(G_t) > 0)){
#   #   scalingFactor <- diag(G_t)^(-.5)  
#   # } else {
#   #   scalingFactor <- c(1,1)
#   # }
#   
#   Glambda <- G_lambda(a_l, b_l)
#   
#   nat_grad <- solve(Glambda) %*% grad_est
#   
#   # a_l_T <-  a_l_T + rho_l_T * grad_est[1]
#   # b_l_T <-  b_l_T + rho_l_T * grad_est[2]
#   a_l <-  a_l + rho * nat_grad[1]
#   b_l <-  b_l + rho * nat_grad[2]
#   # a_l <-  a_l + rho * scalingFactor[1] * grad_est[1]
#   # b_l <-  b_l + rho * scalingFactor[2] * grad_est[2]
#   
#   G_t <- G_t + grad_est %*% t(grad_est)
#   
#   list("a_l" = a_l,
#        "b_l" = b_l,
#        "G_t" = G_t)
# }
# 
# update_sig_T <- function(y, a_l, b_l, 
#                        a_sig, b_sig,
#                        sig2, x_gp,
#                        a_sig_0, b_sig_0,
#                        G_t,
#                        sims_grad, rho){
#   
#   grad_est <- apply(sapply(1:sims_grad, function(s){
#     sig_s <- rgamma(1, a_sig, b_sig)
#     l_s <- rgamma(1, a_l, b_l)
#     Grad_logq(sig_s, a_sig, b_sig) %*% (logp_T(y, l_s, sig_s, x_gp, sig2, a_sig_0, b_sig_0) - 
#                                           logq(sig_s, a_sig, b_sig))
#   }), 1, mean)
#   
#   if(all(diag(G_t) > 0)){
#     scalingFactor <- diag(G_t)^(-.5)  
#   } else {
#     scalingFactor <- c(1,1)
#   }
#   
#   Glambda <- G_lambda(a_sig, b_sig)
#   
#   nat_grad <- solve(Glambda) %*% grad_est
#   
#   a_sig <-  a_sig + rho * nat_grad[1]
#   b_sig <-  b_sig + rho * nat_grad[2]
#   # a_sig <-  a_sig + rho * scalingFactor[1] * grad_est[1]
#   # b_sig <-  b_sig + rho * scalingFactor[2] * grad_est[2]
#   
#   G_t <- G_t + grad_est %*% t(grad_est)
#   
#   list("a_sig" = a_sig,
#        "b_sig" = b_sig,
#        "G_t" = G_t)
# }
# 
# update_l_S <- function(y, a_l, b_l, 
#                      a_sig, b_sig,
#                      sig2, x_gp,
#                      a_l_0, b_l_0,
#                      sims_grad, rho){
#   
#   grad_est <- apply(sapply(1:sims_grad, function(s){
#     l_s <- rgamma(1, a_l, b_l)
#     sig_s <- rgamma(1, a_sig, b_sig)
#     Grad_logq(l_s, a_l, b_l) %*% (logp_s(y, l_s, sig_s, x_gp, sig2) - logq(l_s, a_l, b_l))
#   }), 1, mean)
#   
#   Glambda <- G_lambda(a_l, b_l)
#   
#   nat_grad <- solve(Glambda) %*% grad_est
#   
#   # a_l_T <-  a_l_T + rho_l_T * grad_est[1]
#   # b_l_T <-  b_l_T + rho_l_T * grad_est[2]
#   a_l <-  a_l + rho * nat_grad[1]
#   b_l <-  b_l + rho * nat_grad[2]
#   
#   list("a_l" = a_l,
#        "b_l" = b_l)
# }
# 
# update_sig_S <- function(y, a_l, b_l, 
#                        a_sig, b_sig,
#                        sig2, x_gp,
#                        a_sig_0, b_sig_0,
#                        sims_grad, rho){
#   
#   # grad_est <- apply(sapply(1:sims_grad, function(s){
#   #   sig_s <- rgamma(1, a_sig, b_sig)
#   #   l_s <- rgamma(1, a_l, b_l)
#   #   Grad_logq(sig_s, a_sig, b_sig) %*% (logp_s(y, l_s, sig_s, x_gp, sig2) - 
#   #                                         logq(sig_s, a_sig, b_sig))
#   # }), 1, mean)
#   
#   grad_sims <- sapply(1:sims_grad, function(s){
#     l_s <- rgamma(1, a_l_S, b_l_S)
#     sig_s <- rgamma(1, a_sigma_S, b_sigma_S)
#     h_i <- Grad_logq(sig_s, a_sigma_S, b_sigma_S)
#     d_i <- logp_s(y_true,
#                   l_s, sig_s^2, X, 0) -
#       logq(sig_s, a_sigma_S, b_sigma_S)
#     f_i <- h_i %*% d_i
#     cbind(h_i, rep(d_i, length(f_i)), f_i)
#   })
#   
#   a_i <- sum(sapply(1:2, function(d){
#     cov(grad_sims[d + 4,], grad_sims[d,])
#   })) /  sum(sapply(1:2, function(d){
#     var(grad_sims[d,])
#   }))
#   
#   grad_sims_new <- sapply(1:sims_grad, function(s){
#     grad_sims[1:2, s] * (grad_sims[3,s] - a_i)   
#   })
#   
#   grad_est <- apply(grad_sims_new, 1, mean)
#   
#   Glambda <- G_lambda(a_sig, b_sig)
#   
#   nat_grad <- solve(Glambda) %*% grad_est
#   
#   a_sig <-  a_sig + rho * nat_grad[1]
#   b_sig <-  b_sig + rho * nat_grad[2]
#   
#   list("a_sig" = a_sig,
#        "b_sig" = b_sig)
# }

# ----------------

rinvgamma <- function(a, b){
  1 / stats::rgamma(1, shape = a, rate = b)
}


dinvgamma <- function (x, shape, scale = 1)  {
  log.density <- shape * log(scale) - lgamma(shape) - 
    (shape + 1) * log(x) - (scale /  x)
  return(exp(log.density))
}


update_sigma <- function(a_s, k_s, a_sigma, b_sigma){
  
  a_s_sites <- a_s[!duplicated(k_s$Site)]
  n <- length(a_s_sites)
  
  rinvgamma(a_sigma + (n / 2), b_sigma + sum(a_s_sites^2)/2)
  
}

update_l_sigma_integrated <- function(l_T, sigma_T, a_l_T, b_l_T, a_ig, b_ig, Y, sigma_psi,
                                      z, XbetaY, Xbetas, Xbeta_cov, eps_s, Omega,  
                                      b_psi_b, sd_l, sd_sigma_T,  X_y_index, usingSpatial){
  
  l_T_star <- stats::rnorm(1, l_T, sd_l)
  sigma_T_star <- stats::rnorm(1, sigma_T, sd_sigma_T)
  
  if(l_T_star < 4 & l_T_star > 0 & sigma_T_star > 0){
    
    K_l <- K(1:Y, 1:Y, sigma_T^2, l_T) + sigma_psi^2 + diag(exp(-10), nrow = Y)
    K_l_star <- K(1:Y, 1:Y, sigma_T_star^2, l_T_star) + sigma_psi^2 + diag(exp(-10), nrow = Y)
    
    # term 1
    
    tXOmegX <- matrixProductXtOmegaX_year(Y, Omega, X_y_index)
    
    logterm1_1 <- -.5 * log(FastGP::rcppeigen_get_det(K_l_star %*% tXOmegX + diag(1, nrow = Y)))
    logterm1_2 <- -.5 * log(FastGP::rcppeigen_get_det(K_l %*% tXOmegX + diag(1, nrow = Y)))
    
    logdeterminantsRatio <- logterm1_1 - logterm1_2
    
    # term 2
    
    logmuKlmu <- - .5 * b_psi_b %*% (FastGP::rcppeigen_invert_matrix(K_l_star) - 
                                       FastGP::rcppeigen_invert_matrix(K_l)) %*% b_psi_b
    
    # term 3
    
    if(usingSpatial){
      c_i <- Xbetas + Xbeta_cov + eps_s  
    } else {
      c_i <- Xbeta_cov + eps_s
    }
    
    
    k_i <- z - .5
    z_tilde <- k_i - c_i * as.vector(Omega)
    
    Xtz_tilde <- XpsiYz(X_y_index, z_tilde, Y)
    
    mu_tilde <- Xtz_tilde + FastGP::rcppeigen_invert_matrix(K_l) %*% b_psi_b
    mu_tilde_star <- Xtz_tilde + FastGP::rcppeigen_invert_matrix(K_l_star) %*% b_psi_b
    
    XtOXpB <- tXOmegX + FastGP::rcppeigen_invert_matrix(K_l)
    XtOXpB_star <- tXOmegX + FastGP::rcppeigen_invert_matrix(K_l_star)
    
    log_term2ratio_star <- .5 * ( t(mu_tilde_star) %*% FastGP::rcppeigen_invert_matrix(XtOXpB_star) %*% mu_tilde_star )
    log_term2ratio <- .5 * ( t(mu_tilde) %*% FastGP::rcppeigen_invert_matrix(XtOXpB) %*% mu_tilde )
    
    logterm3 <- log_term2ratio_star - log_term2ratio
    
    likelihood <- exp(logdeterminantsRatio + logmuKlmu +  logterm3)
    
    prior_l <- exp(stats::dgamma(l_T_star, a_l_T, b_l_T, log = T) - stats::dgamma(l_T, a_l_T, b_l_T, log = T))
    prior_sigma <- exp(dinvgamma(sigma_T_star, a_ig, b_ig) - dinvgamma(sigma_T, a_ig, b_ig))
    
    posterior <-  likelihood * prior_l * prior_sigma
    
    if(!is.na(posterior)){
      if(stats::runif(1) < posterior){
        l_T <- l_T_star
        sigma_T <- sigma_T_star
      }  
    }
    
  }
  
  list("l_T" = l_T, "sigma_T" = sigma_T)
}

rcpp_log_dmvnorm_fast <- function (inv_S, diag_S, sigma_s, x) {
  n <- length(x)
  # return((-n/2) * log(2 * pi) - sum(log(diag_S)) - (1/2) * x %*% inv_S %*% x)
  return((-n/2) * log(2 * pi) - sum(log(diag_S) + log(sigma_s)) -
           (1/2) * (1 / sigma_s^2) * x %*% inv_S %*% x)
}

sample_l_grid <- function(l_s_grid, inv_K_s_grid, diag_K_s_grid,
                          a_l_s, b_l_s, a_s, sigma_s){
  
  posterior_val <- rep(NA, length(l_s_grid))
  
  for (j in 1:length(l_s_grid)) {
    
    l_s <- l_s_grid[j]
    # K_s_grid_j <- K_s_grid[,,j] * sigma_s^2
    # inv_K_s_grid_j <- inv_K_s_grid[,,j] / sigma_s^2
    # diag_K_s_grid_j <- diag_K_s_grid[,j] * sigma_s
    
    # loglikelihood <- rcpp_log_dmvnorm_fast(inv_K_s_grid[,,j],
    #                                        diag_K_s_grid[,j], sigma_s, a_s)
    
    loglikelihood <- ((-length(a_s)/2) * log(2 * pi) - sum(log(diag_K_s_grid[,j]) + log(sigma_s)) -
                        (1/2) * (1 / sigma_s^2) * a_s %*% inv_K_s_grid[,,j] %*% a_s)
    
    # loglikelihood <- rcpp_log_dmvnorm_fast(1, inv_K_s_grid_j,
    # diag_K_s_grid_j,  a_s)
    
    # Sigma_l <- K2(X_tilde, X_tilde, sigma_s^2, l_s) + diag(exp(-10), nrow = nrow(X_tilde))
    
    # (loglikelihood2 <- rcpp_log_dmvnorm( Sigma_l, rep(0, X_centers), a_s, F))
    
    logPrior <- dgamma(l_s, a_l_s, b_l_s, log = T)
    
    posterior_val[j] <- logPrior + loglikelihood
    
  }
  
  # posterior_val <- sample_l_grid_cpp(l_s_grid, sigma_s, inv_K_s_grid, diag_K_s_grid,
  # a_l_s, b_l_s, a_s)
  
  posterior_val <- posterior_val - max(posterior_val)
  
  index_l_grid <- sample(1:length(l_s_grid), size = 1, 
                         prob = exp(posterior_val) / sum(exp(posterior_val)))
  
  index_l_grid
}

update_hyperparameters <- function(l_T, a_l_T, b_l_T, sd_l_T, sd_sigma_T,
                                   sigma_T, a_sigma_T, b_sigma_T, Y,
                                   beta_psi, inv_B_psi, 
                                   b_psi, sigma_psi,
                                   l_s_grid, K_s_grid, inv_K_s_grid, 
                                   inv_chol_K_s_grid, diag_K_s_grid,
                                   a_l_s, b_l_s, 
                                   sigma_s, a_sigma_s, b_sigma_s, X_tilde,
                                   a_sigma_eps, b_sigma_eps,
                                   usingSpatial, 
                                   XbetaY, Xbetas, Xbeta_cov,
                                   eps_s, k_s, 
                                   z, X_psi, Omega, X_y_index){
  
  
  # update l_T and sigma_T -------------------
  
  list_l_sigma <- update_l_sigma_integrated(l_T, sigma_T, a_l_T, b_l_T, a_sigma_T, b_sigma_T, Y, 
                                            sigma_psi,
                                            z, betaY, Xbetas, Xbeta_cov, eps_s, Omega,  
                                            b_psi[1:Y], sd_l = sd_l_T, 
                                            sd_sigma_T, X_y_index, usingSpatial)
  l_T <- list_l_sigma$l_T
  sigma_T <- list_l_sigma$sigma_T
  
  # reupdate inv_B_psi ------------
  
  K_l <- K(1:Y, 1:Y, sigma_T^2, l_T) + sigma_psi^2 + diag(exp(-8), nrow = Y)
  
  inverse_Kl <- FastGP::rcppeigen_invert_matrix(K_l)
  inverse_Kl[lower.tri(K_l)] <- t(inverse_Kl)[lower.tri(K_l)]
  
  inv_B_psi[1:Y, 1:Y] <- inverse_Kl
  
  # update l_s -------------------
  
  if(usingSpatial){
    index_ls_grid <- sample_l_grid(l_s_grid, inv_K_s_grid, diag_K_s_grid,
                                   a_l_s, b_l_s, beta_psi[Y + 1:nrow(X_tilde)], sigma_s)
    l_s <- l_s_grid[index_ls_grid]
  } else {
    l_s <- 0
    index_ls_grid <- 0
  }
  
  # update sigma_s ---------------------
  
  if(usingSpatial){
    
    X_centers <- nrow(X_tilde)
    
    inv_chol_Kl <- inv_chol_K_s_grid[,,index_ls_grid]
    
    a_s <- beta_psi[Y + 1:X_centers]
    
    a_ig <- X_centers / 2
    Ltas <- inv_chol_Kl %*% a_s
    b_ig <- (t(Ltas) %*% Ltas) / 2
    
    sigma_s <- sqrt(rinvgamma(a_sigma_s + a_ig, b_sigma_s + b_ig))
  } 
  
  # reupdate inv_B_psi ------------
  
  if(usingSpatial){
    X_centers <- nrow(X_tilde)
    
    inv_K_lsigma <- inv_K_s_grid[,,index_ls_grid] / sigma_s^2
    
    inv_B_psi[Y + 1:X_centers, Y + 1:X_centers] <- inv_K_lsigma  
  }
  
  # update sigma_as ------------------------------
  
  sigma_eps <- sqrt(update_sigma(eps_s, k_s, a_sigma_eps, b_sigma_eps))
  
  list("l_T" = l_T,
       "sigma_T" = sigma_T,
       "inv_B_psi" = inv_B_psi,
       "l_s" = l_s,
       "index_ls_grid" = index_ls_grid,
       "sigma_s" = sigma_s,
       "sigma_eps" = sigma_eps)
}
