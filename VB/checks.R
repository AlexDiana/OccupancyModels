# YEAR EFFECT ---------

years <- modelResults$dataCharacteristics$Years
Y <- length(years)

mu_beta_psi <- modelResults$modelOutput$mu_beta_psi
Sigma_beta_psi <- modelResults$modelOutput$Sigma_beta_psi

CI_yearseffect <- sapply(1:Y, function(j){
  logit(mu_beta_psi[j] + c(-1,1,0) * 1.96 * sqrt(Sigma_beta_psi[j,j]))
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

# COVARIATES OCCUPANCY ----------

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

# beta_psi_output <- modelResults$modelOutput$beta_psi_output

# beta_psi_output <- apply(beta_psi_output, 3, c)

betaeffect <- sapply(1:(ncov_psi + numTimeSpaceCov), function(k){
  mu_beta_psi[Y + X_centers + k] + c(-1,1,0) * 1.96 * sqrt(Sigma_beta_psi[Y + X_centers + k,
                                                                          Y + X_centers + k])
}) 
#   apply(beta_psi_output[,Y + X_centers + 1:ncov_psi_all,drop = F], 2, function(x){
#   quantile(x, probs = c(0.025,0.5,0.975))
# })
# betaeffect <- apply(beta_psi_output[,Y + X_centers + 1:ncov_psi_all,drop = F], 2, function(x){
#   quantile(x, probs = c(0.025,0.5,0.975))
# })

beta_psi_true_all <- beta_psi_true
if(usingSpatial){
  beta_psi_true_all <- c(beta_psi_tsp_true, beta_psi_true_all)
}

betaeffect <- rbind(betaeffect, beta_psi_true_all)
colnames(betaeffect) <- namesCovariates_all

ggplot(data = NULL, aes(x = colnames(betaeffect),
                        y = betaeffect[3,],
                        ymin = betaeffect[1,],
                        ymax = betaeffect[2,])) + geom_point() + geom_errorbar() +
  geom_point(data = NULL, aes(x = colnames(betaeffect),
                              y = betaeffect[4,]), color = "red", size = 4) + 
  ylab("Effect") + scale_x_discrete(name = "Covariates") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 16, face = "bold"),
        panel.grid.major = element_line(colour="grey", size=0.015),
        panel.background = element_rect(fill = "white", color = "black"))

# DETECTION PROB --------

mu_beta_p <- modelResults$modelOutput$mu_beta_p
Sigma_beta_p <- modelResults$modelOutput$Sigma_beta_p

if(usingYearDetProb){
  CI_yearseffect <- sapply(1:Y, function(j){
    logit(mu_beta_p[j] + c(-1,1,0) * 1.96 * sqrt(Sigma_beta_p[j,j]))
  })
} else {
  CI_yearseffect <- logit(mu_beta_p[1] + c(-1,1,0) * 1.96 * sqrt(Sigma_beta_p[1,1]))
}


prob_true <- logit(mu_p_true #+ 
                   # beta_psi_t_true * sort(unique(modelResults$dataCharacteristics$X_psi_yearcov)) +
                   # mean(a_s_site_true)
)

if(usingYearDetProb){
  
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
} else {
  ggplot(data = NULL, aes(x = 1,
                          y = CI_yearseffect[3],
                          ymin = CI_yearseffect[1],
                          ymax = CI_yearseffect[2])) + 
    geom_point() + geom_errorbar() + geom_line() + 
    geom_point(data = NULL, aes(x = 1, y = prob_true), size = 2, color = "red") +
    ylab("Occupancy probability") + scale_x_continuous(name = "Year", breaks = years) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black")) + 
    ylim(c(.1,.3))
}

# COVARIATES DETECTION ----------

X_tilde <- modelResults$dataCharacteristics$X_tilde
if(usingSpatial){
  X_centers <- nrow(X_tilde)  
} else {
  X_centers <- 0
}

years <- modelResults$dataCharacteristics$Years
Y <- length(years)

p_intercepts <- ifelse(usingYearDetProb, Y, 1)

numTimeSpaceCov <- modelResults$dataCharacteristics$numTimeSpaceCov

namesCovariates <- modelResults$dataCharacteristics$nameVariables_p
ncov_p <- length(namesCovariates)

# beta_psi_output <- modelResults$modelOutput$beta_psi_output

# beta_psi_output <- apply(beta_psi_output, 3, c)

betaeffect <- sapply(1:ncov_p, function(k){
  mu_beta_p[p_intercepts + k] + c(-1,1,0) * 1.96 * sqrt(Sigma_beta_p[p_intercepts + k,
                                                                     p_intercepts + k])
}) 
#   apply(beta_psi_output[,Y + X_centers + 1:ncov_psi_all,drop = F], 2, function(x){
#   quantile(x, probs = c(0.025,0.5,0.975))
# })
# betaeffect <- apply(beta_psi_output[,Y + X_centers + 1:ncov_psi_all,drop = F], 2, function(x){
#   quantile(x, probs = c(0.025,0.5,0.975))
# })

beta_psi_true_all <- beta_psi_true
if(usingSpatial){
  beta_psi_true_all <- c(beta_psi_tsp_true, beta_psi_true_all)
}

betaeffect <- rbind(betaeffect, beta_p_true)
colnames(betaeffect) <- namesCovariates

ggplot(data = NULL, aes(x = colnames(betaeffect),
                        y = betaeffect[3,],
                        ymin = betaeffect[1,],
                        ymax = betaeffect[2,])) + geom_point() + geom_errorbar() +
  geom_point(data = NULL, aes(x = colnames(betaeffect),
                              y = betaeffect[4,]), color = "red", size = 4) + 
  ylab("Effect") + scale_x_discrete(name = "Covariates") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 16, face = "bold"),
        panel.grid.major = element_line(colour="grey", size=0.015),
        panel.background = element_rect(fill = "white", color = "black"))

# PLOT OF L_T ----------

l_T_output <- modelResults$modelOutput$l_T

ggplot(data = NULL) + stat_function(fun = dgamma, args = list(shape = a_l_T,
                                                              rate = b_l_T)) + 
  xlim(c(0, 10)) + theme_bw()

a_l_T / b_l_T

# PLOT OF SIGMA_T -------------

a_sigma_T <- modelResults$modelOutput$a_sigma_T
b_sigma_T <- modelResults$modelOutput$b_sigma_T

ggplot(data = NULL) + stat_function(fun = dgamma, args = list(shape = a_sigma_T,
                                                              rate = b_sigma_T)) + 
  xlim(c(0, 10)) + theme_bw()

a_sigma_T / b_sigma_T

# SPATIAL ------


X_tilde <- modelResults$dataCharacteristics$X_tilde
X_centers <- nrow(X_tilde)

years <- modelResults$dataCharacteristics$Years
Y <- length(years)


mu_beta_psi <- modelResults$modelOutput$mu_beta_psi

# eps_unique_output <- apply(eps_unique_output, 2, c)
# beta_psi_output <- apply(beta_psi_output, 3, c)

a_s_star <- mu_beta_psi[Y + 1:nrow(X_tilde)]

# sapply(1:nrow(X_tilde), function(i){
#   mean( )#+ eps_unique_output[i])
#   # )
# }) 

ggplot() + geom_point(data = NULL, aes(x = X_tilde[,1],
                                       y = X_tilde[,2], color = a_s_star), 
                      size = 8, shape = 15) +
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
