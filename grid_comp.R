library(mvtnorm)

l_s <- .1
sigma_s <- 1

eps <- exp(-5)

n_1 <- 40^2
n_2 <- 20^2
X_1 <- expand.grid(seq(0, 1, length.out = sqrt(n_1)),
                 seq(0, 1, length.out = sqrt(n_1)))
X_1 <- scale(X_1)
X_2 <- expand.grid(seq(0, 1, length.out = sqrt(n_2)),
                 seq(0, 1, length.out = sqrt(n_2)))
X_2 <- scale(X_2)

nuggetMat <- diag(exp(-10), nrow = n_2)

X1 <- K2(X_1, X_2, sigma_s^2, l_s)
X2 <- K2(X_2, X_2, sigma_s^2, l_s)
X1_star <- K2(X_1, X_1, sigma_s^2, l_s)
# X <- X1 %*% solve(X2 + nuggetMat) %*% t(X1)

# a_tilde <- rcpp_rmvnorm(1, X2, rep(0, n_2))#mvrnorm(1, rep(0, n_2), X2)
# a_s <- as.vector(X1 %*% solve(X2 + nuggetMat) %*% t(a_tilde))
a_s <- as.vector(rcpp_rmvnorm(1, X1_star, rep(0, n_1)))#mvrnorm(1, rep(0, n_1), X1_star)

# plot 
{
  ggplot() +
    geom_point(data = NULL, aes(x = X_1[,1],
                                y = X_1[,2], alpha = a_s), size = 3, shape = 15) +
    xlab("X") + ylab("Y") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          # legend.position = "none",
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          # panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black"))
    
}

log_dgen <- function(x, sigma_s, l_s){
  
  X1 <- K2(X_1, X_2, 1^2, l_s)
  X2 <- K2(X_2, X_2, 1^2, l_s)
  
  list_params <- ginvsquare_pseudodet(X1, X2, eps)
  sq_current <- list_params$sq
  lpdet <- list_params$lpdet
  nval <- list_params$nval
  
  xsq <- x %*% sq_current
  
  # loglikelihood <- - nval / 2 * (log(2 * pi) + log(sigma_s^2)) - .5 * ldet_grid[j] -
  - nval / 2 * (log(2 * pi) + log(sigma_s^2)) - .5 * lpdet -
    (1/2) *  (1/ sigma_s^2) * (xsq %*% t(xsq))
  # - nval / 2 * (log(2 * pi) ) - .5 * lpdet -
  #   (1/2) * (xsq %*% t(xsq))
  
}

log_dgen2 <- function(x, sigma_s, l_s){
  
  X1 <- K2(X_1, X_2, sigma_s, l_s)
  X2 <- K2(X_2, X_2, sigma_s, l_s)
  
  det(X2)
  
  X <- X1 %*% solve(X2 + diag(exp(-25), nrow = n_2))
  
  a_s_tilde <- lm(a_s ~ X - 1)$coeff
  max(a_s_tilde)
  # log(det(K2(X_2, X_2, sigma_s^2, l_s) +
  #           diag(eps, nrow = n_2)))
  
  FastGP::rcpp_log_dmvnorm(S = (K2(X_2, X_2, sigma_s^2, l_s) + 
                                  diag(eps, nrow = n_2)), mu = rep(0, n_2),
                           a_s_tilde, istoep = F)
  
  # a_s_tilde %*% rcppeigen_invert_matrix(K2(X_2, X_2, sigma_s^2, l_s) + 
  #                           diag(eps, nrow = n_2)) %*% a_s_tilde
  
}

log_dmvnorm <- function(x, sigma_s, l_s){
  
  FastGP::rcpp_log_dmvnorm(S = (K2(X_1, X_1, sigma_s^2, l_s) + 
                                  diag(eps, nrow = n_1)), mu = rep(0, n_1),
                          x, istoep = F)
  
}

# sigma 

sigma_s_grid <- seq(sigma_s / 1.2, sigma_s * 1.1, length.out = 10)
i <- 1
sigma_s_values <- sapply(sigma_s_grid, function(x){
  print(i)
  i <<- i + 1
  # log_dmvnorm(a_s, x, l_s)
  # log_dgen(a_s, x, l_s)
  log_dgen2(a_s, x, l_s)
})

sigma_s_values <- sigma_s_values - max(sigma_s_values)
sigma_s_values <- exp(sigma_s_values) / sum(exp(sigma_s_values))

qplot(sigma_s_grid, sigma_s_values) + 
  geom_vline(aes(xintercept = sigma_s))

# l 

l_s_grid <- seq(l_s / 2, l_s * 1.2, length.out = 10)
i <- 1
l_s_values <- sapply(l_s_grid, function(x){
  print(i)
  i <<- i + 1
  # log_dmvnorm(a_s, x, l_s)
  log_dgen(a_s, sigma_s, x)
  # log_dgen2(a_s, sigma_s, x)
})

l_s_values <- l_s_values - max(l_s_values)
l_s_values <- exp(l_s_values) / sum(exp(l_s_values))

qplot(l_s_grid, l_s_values) + 
  geom_vline(aes(xintercept = l_s))
