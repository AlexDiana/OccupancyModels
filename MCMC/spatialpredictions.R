
X_tilde <- modelResults_MCMC$dataCharacteristics$X_tilde
X_centers <- nrow(X_tilde)
years <- modelResults_MCMC$dataCharacteristics$Years
Y <- length(years)

beta_psi_output <- modelResults_MCMC$modelOutput$beta_psi_output
beta_psi_output <- apply(beta_psi_output, 3, c)

l_s_output <- modelResults_MCMC$modelOutput$l_s_output
l_s_mean <- mean(l_s_output)

x_grid <- seq(min(X_tilde[,1]), max(X_tilde[,1]), length.out = 60)
y_grid <- seq(min(X_tilde[,2]), max(X_tilde[,2]), length.out = 60)
siteLocations <- as.matrix(expand.grid(x_grid, y_grid))

# prediction at new points
{
  K_staru <- K2(siteLocations, X_tilde, 1, l_s_mean)
  inv_K_uu <- solve(K2(X_tilde, X_tilde, 1, l_s_mean) + 
                      diag(exp(-10), nrow = nrow(X_tilde)))
  X_coeff <- K_staru %*% inv_K_uu
  
  # K_biggest <- apply(K_staruKu_all[,,i], 1, function(x){
  #   -sort(-x)[1:maxPoints]
  # })
  # 
  # which_K_biggest <- apply(K_staruKu_all[,,i], 1, function(x){
  #   order(-x)[1:maxPoints]
  # })
  # 
  # X_s_index_all[,,i] <- t(which_K_biggest)
  # for (j in 1:nrow(X_psi)) {
  #   K_staruKu_all[j,-X_s_index_all[j,,i],i] <- 0
  #   K_staruKu_all[j,,i] <- K_staruKu_all[j,,i] / sum(K_staruKu_all[j,,i])
  # }
  
  beta_mean <- apply(beta_psi_output, 2, mean)
  
  a_s_star <- X_coeff %*% beta_mean[Y + seq_len(X_centers)]
}

npoints <- nrow(siteLocations)
datapoly <- data.frame(
  id = rep(1:npoints, each = 4),
  x = rep(siteLocations[,1], each =  4),
  y = rep(siteLocations[,2], each =  4),
  value = rep(a_s_star, each =  4)
)

# newGridStep <- (unique(sort(siteLocation[,1]))[2] - 
# unique(sort(siteLocation[,1]))[1] ) / 1.2
newGridStep <- unique(sort(siteLocations[,1]))[2] - unique(sort(siteLocations[,1]))[1] 

variations_x <- rep(c(newGridStep, -newGridStep,
                      -newGridStep, newGridStep), times = npoints)
variations_y <- rep(c(-newGridStep, -newGridStep,
                      newGridStep, newGridStep), times = npoints)

datapoly$x <- datapoly$x + variations_x
datapoly$y <- datapoly$y + variations_y

ggplot(datapoly, aes(x = x, y = y)) +
  geom_polygon(aes(alpha = value, group = id), fill = "black") +
  scale_alpha_continuous(range = c(0, 0.8)) +
  xlab("") + ylab("") +
  scale_x_continuous(breaks = c()) +
  scale_y_continuous(breaks = c()) +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "none",
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 13, face = "bold", angle = 90),
        # panel.grid.major = element_line(colour="grey", size=0.015),
        panel.background = element_rect(fill = "white", color = "black"))

