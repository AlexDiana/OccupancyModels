ELBO <- function(y, 
                 mu_beta_psi, Sigma_beta_psi,
                 mu_beta_p, Sigma_beta_p,
                 gamma_z,
                 l_T, 
                 a_sigma_T, b_sigma_T,
                 usingSpatial,
                 l_S, 
                 a_sigma_S, b_sigma_S){
  
  nsims <- 40
  elbovals <- sapply(1:nsims, function(i){
    
    # parameter simulations
    {
      z <- simulate_q_z(gamma_z)
      z_all <- z[indexes_occobs]
      # z <- sapply(1:length(gamma_z), function(j){
      #   rbinom(1, 1, gamma_z[j])
      # })
      # z_all <- z[indexes_occobs]
      
      sigmasq_T <- rinvgamma(a_sigma_T, b_sigma_T)
      
      if(usingSpatial){
        
        sigmasq_S <- rinvgamma(a_sigma_S, b_sigma_S)
        
      }
      
      beta_psi <- mvrnormArma(mu_beta_psi, Sigma_beta_psi)
      eps_s <- simulate_eps(mu_eps, sd_eps)
      
      # eps_s <- sapply(1:length(mu_eps), function(j){
      #   rnorm(1, mu_eps[j], sd_eps[j])
      # })
      Xpsibeta <- tXmubetapsi(X_psi, X_y_index, X_s_index,
                                   beta_psi, Y, X_centers,
                                        ncov_psi + numTimeSpaceCov) + eps_s
      # Xpsibeta <- beta_psi[X_y_index] + beta_psi[Y + X_s_index] + 
        # X_psi[,Y + X_centers + 
                # 1:(numTimeSpaceCov + ncov_psi)] %*% beta_psi[Y + X_centers + 
                                                               # 1:(numTimeSpaceCov + ncov_psi)] + eps_s
      psi <- logit(Xpsibeta)
      
      beta_p <- mvrnorm(1, mu_beta_p, Sigma_beta_p)
      
      p <- logit(X_p %*% beta_p)  
    }
    
    # log p(x, z)
    {
      # y | z, p
      
      y_likelihood <- y_lik(y, p, z_all)
      
      # z | psi
      z_lik <- log_q_z(psi, z)
        
      # beta_psi | l_T, sigma_T
      Sigma_psi_T <- K(1:Y, 1:Y, sigma_T, l_T)
      (lik_betapsi_T <- log_dmvnorm_cpp(beta_psi[1:Y], rep(0, Y), Sigma_psi_T))
      
      # beta_psi | l_S, sigma_S
      if(usingSpatial){
        Sigma_psi_S <- K2(X_tilde, X_tilde, sigma_S, l_S)
        lik_betapsi_S <- rcpp_log_dmvnorm(Sigma_psi_S, 
                                          rep(0, X_centers), 
                                          beta_psi[Y + 1:X_centers], F)  
      } else {
        lik_betapsi_S <- 0
      }
      
      # lik_betapsi_S <- log_dmvnorm_cpp(beta_psi[Y + 1:X_centers], 
      #                                  rep(0, X_centers), 
      #                                  Sigma_psi_S)
      

      logpxz <- y_likelihood + z_lik + lik_betapsi_T + lik_betapsi_S
      
      
      c(y_likelihood, z_lik, lik_betapsi_T, lik_betapsi_S)
    }
    
    # log q
    {
     
      sum_log_q_z <- log_q_z(gamma_z, z)
        
        # sapply(1:length(gamma_z), function(j){
        # dbinom(z[j], 1, gamma_z[j])
      # })
      
      log_beta_psi <- log_dmvnorm_cpp(beta_psi, mu_beta_psi, Sigma_beta_psi)
      
      sum_log_eps_s <- log_q_eps(eps_s, mu_eps, sd_eps)
      
      # log_eps_s <- sapply(1:length(mu_eps), function(j){
      #   dnorm(eps_s[j], mu_eps[j], sd_eps[j], log = T)
      # })
      
      log_beta_p <- log_dmvnorm_cpp(beta_p, mu_beta_p, Sigma_beta_p)
      
      logq <- sum_log_q_z + log_beta_psi + sum_log_eps_s + 
        log_beta_p
    }
    
    logpxz - logq
    
  })
  
  mean(elbovals)
  
}

#' Model fitting
#' 
#' @description This function fits the occupancy model of Diana et al. (2021). 
#' Note that in the following the parameters are described with the notations used in the paper.
#' 
#' @param data The data frame containing the data.
#' @param index_year The index of the column containing the year variable.
#' @param index_site The index of the column containing the site variable.
#' @param index_occ The index of the column containing the detections.
#' @param index_spatial_x The index of the x coordinate of the site.
#' @param index_spatial_y The index of the y coordinate of the site. 
#' @param covariates_psi_text Indexes of the column of the occupancy probability. To be separated by a comma, eg. "5,6,8". Set to "0" if no covariate is available
#' @param covariates_p_text Indexes of the column of the detection probability. 
#' @param prior_psi Prior mean for the occupancy probability
#' @param sigma_psi Standard deviation for the prior on the occupancy probability
#' @param prior_p Prior mean for the detection probability
#' @param sigma_p Standard deviation for the prior on the detection probability
#' @param usingYearDetProb Should the model include year-specific detection probabilities (as opposed to a constant one)?
#' @param usingSpatial Should the model include the auto-correlated spatial effects?
#' @param gridStep Step of the grid to use for the approximation of the auto-correlated spatial effects. Use \code{\link{buildSpatialGrid}} to show the grid for a value of \code{gridStep}.
#' @param storeRE Should the model store the site-specific independent random effects for each iteration (instead of just their mean across all chain)? Not suggested if the number of sites is greater than 1000.
#' @param nchain Number of chains.
#' @param nburn Number of burn-in iterations.
#' @param niter Number of (non burn-in) iterations.
#' @param verbose Should the progress of the MCMC be printed?.
#' @param computeGOF Should the model perform calculations of the goodness of fit?
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return 
#' A list with components:
#' 
#' \itemize{
#' 
#'  \item \code{modelResults} A list with components:
#'       
#'       \describe{
#'       
#'       \item{\code{beta_psi_output}}{ An array of dimensions \code{nchain} x \code{niter} x (number of coefficients for 
#'       occupancy probability) containing the values of the coefficients of the occupancy probability.
#'      The coefficients are reported in the order (year r.e., space r.e., covariates for time-space interaction,
#'       standard covariates).}
#'       
#'       \item{\code{eps_unique_output}}{ if \code{storeRE = T}, an array of dimensions 
#'       \code{nchain} x \code{niter} x S containing the values of the site-specific independent random effects
#'        of the occupancy probability. If  \code{storeRE = F}, a matrix of dimensions \code{nchain} x S 
#'        containing the mean value of the random effect across a chain.}
#'       
#'       \item{\code{beta_p_output}}{ An array of dimensions \code{nchain} x \code{niter} x (number of coefficients for 
#'       detection probability) containing the values of the coefficients of the detection probability.
#'        The coefficients are reported in the order (intercepts, covariates).}
#'       
#'        \item{\code{sigma_T_output}}{ A array of dimensions \code{nchain} x \code{niter} containing the values of
#'        \eqn{\sigma_T}.}
#'        
#'        \item{\code{l_T_output}}{ A array of dimensions \code{nchain} x \code{niter} containing the values of
#'        \eqn{l_T}.}
#'        
#'        \item{\code{sigma_s_output}}{ A array of dimensions \code{nchain} x \code{niter} containing the values of
#'        \eqn{\sigma_S}.}
#'        
#'        \item{\code{l_s_output}}{ A array of dimensions \code{nchain} x \code{niter} containing the values of
#'        \eqn{l_S}.}
#'        
#'        \item{\code{sigma_eps_output}}{ A array of dimensions \code{nchain} x \code{niter} containing the values of
#'        \eqn{\sigma_\epsilon}.}
#'        
#'        \item{\code{psi_mean_output}}{ A array of dimensions \code{nchain} x \code{niter} containing the values of
#'        the occupancy index.}
#'        
#'        \item{\code{GOF_output}}{ A list with components:
#'        
#'        \describe{
#'        
#'        \item{  \code{gofYear_output} }{An array of dimensions \code{nchain} x \code{niter} x Y containing the values 
#'          of the test statistics of yearly detections.}
#'        
#'        \item{  \code{gofSpace_output} }{An array of dimensions \code{nchain} x \code{niter} x (M), where M is the number
#'        of regions in the approximation, containing the values of the test statistics of the detections in
#'        each region.}
#'        
#'        \item{  \code{trueYearStatistics} }{ A vector containing the true values of the test statistics of the detections
#'         in each year.}
#'         
#'        \item{  \code{trueSpaceStatistics} }{ A vector containing the true values of the test statistics of the detections
#'         in each region}
#'        
#'        }
#'      
#'        }
#'      
#'      }
#'      
#'        
#'       \item \code{dataCharacteristics}: A list with components:
#'        
#'       \describe{
#'       
#'       \item{\code{Years}}{A vector with the years.}
#'       
#'       \item{\code{X_tilde}}{A matrix of dimension M x 2, where M is the number of points chosen in the 
#'       spatial approximation, with the location of the points used for the approximation. The points are
#'       arranged in the same as order as in the coefficients vector \code{beta_psi_output}.}
#'       
#'       \item{\code{gridStep}}{The width of the grid chosen in the approximation.}
#'       
#'       \item{\code{usingSptial}}{Same as in the input.}
#'       
#'        \item{\code{usingYearDetProb}}{Same as in the input.}
#'        
#'       
#'       }
#'  }
#' 
#' @examples
#' 
#' modelResults <- runModel(sampleData, 
#'                          index_year = 1, 
#'                          index_site = 2, 
#'                          index_occ = 8, 
#'                          index_spatial_x = 3, 
#'                          index_spatial_y = 4, 
#'                          covariates_psi_text = "5", 
#'                          covariates_p_text = "6-7", 
#'                          usingSpatial = TRUE,
#'                          gridStep = .2, 
#'                          nchain = 1, 
#'                          nburn = 100,
#'                          niter = 100)  
#' 
runModel <- function(data, index_year, index_site, 
                     index_occ, index_spatial_x = 0, index_spatial_y = 0,
                     covariates_psi_text, covariates_p_text, 
                     prior_psi = .5, sigma_psi = 2,
                     prior_p = .5, sigma_p = 2, 
                     usingYearDetProb = F,
                     usingSpatial = F, gridStep,
                     storeRE = F, nruns, tol, maxit, 
                     verbose = T, computeGOF = F){
  
  print("Analyzing the data..")
  
  # CLEAN DATA 
  {
    
    colnames(data)[c(index_year, index_site, index_occ)] <- c("Year","Site","Occ") 
    
    if(usingSpatial){
      colnames(data)[c(index_spatial_x, index_spatial_y)] <- c("X_sp","Y_sp")  
    } 
    
    data$Year <- as.numeric(data$Year)
    
    S <- length(unique(data$Site))
    Y <- length(unique(data$Year))
    years <- sort(unique(data$Year))
    
    # transform site codes to numbers
    originalsites <- unique(data$Site)
    names(originalsites) <- seq_along(originalsites)
    data$Site <- as.numeric(names(originalsites)[match(data$Site, originalsites)])
    
    data <- data %>% dplyr::arrange(Site, Year)
    
    # define data for occupancies
    
    {
      k_s <- data %>% dplyr::group_by(Year, Site) %>% 
        dplyr::summarise(Occupied = as.numeric(sum(Occ) > 0),
                         Visits = dplyr::n()) %>% 
        dplyr::arrange(Site, Year) 
      k_s <- as.data.frame(k_s)
      
      # match row in data_p to row in k_s
      indexes_occobs <- rep(1:nrow(k_s), k_s$Visits)
      
      Occs <- data$Occ
      y <- Occs
    }
    
    # covariate for psi
    
    {
      
      X_psi <- data[!duplicated(data[,c("Site","Year")]),] %>% dplyr::arrange(Site, Year)
      
      # year covariates
      {
        X_psi$Year <- as.factor(X_psi$Year)
        X_psi_year <- stats::model.matrix( ~ . - 1, data = X_psi[,c("Year"),drop = F])
        
        X_y_index <- apply(X_psi_year, 1, function(x) {which(x != 0)})
        X_y_index <- unlist(X_y_index)
        
      }
      
      # spatial covariates
      {
        
        if(usingSpatial){
          
          X_sp <- X_psi$X_sp
          Y_sp <- X_psi$Y_sp
          
          uniqueIndexesSite <- which(!duplicated(X_psi$Site))
          X_sp_unique <- X_sp[uniqueIndexesSite]
          Y_sp_unique <- Y_sp[uniqueIndexesSite]
          XY_sp_unique <- cbind(X_sp_unique, Y_sp_unique)
          
          # build the grid
          X_tilde <- as.matrix(buildGrid(XY_sp_unique, gridStep))
          XY_centers <- findClosestPoint(XY_sp_unique, X_tilde)
          
          # get rid of unused cells by refinding the grid
          X_tilde <- X_tilde[sort(unique(XY_centers)),]
          
          # refit
          XY_centers <- findClosestPoint(cbind(X_psi$X_sp,X_psi$Y_sp), X_tilde)
          XY_centers_unique <- findClosestPoint(XY_sp_unique, X_tilde)
          
          X_psi$Center <- as.factor(XY_centers)
          X_s_index <- XY_centers
          
          # X_psi_s <- stats::model.matrix( ~ . - 1, data = X_psi[,c("Center"),drop = F])
          
          X_centers <- nrow(X_tilde)
          
        } else {
          X_centers <- 0
          XY_centers <- NULL
          # X_psi_s <- NULL
          XY_sp_unique <- NULL
          X_s_index <- rep(0, nrow(k_s))
          X_tilde <- NULL
          gridStep <- NULL
          XY_centers_unique <- NULL
        }
        
      }
      
      # time space covariates 
      {
        numTimeSpaceCov <- 0
        
        X_psi_yearcov <- scale(as.numeric(X_psi$Year))
        X_psi_yearcov_values <- sort(unique(X_psi_yearcov))
        
        if(usingSpatial){
          
          X_psi_spcov <- cbind(X_sp, Y_sp)
          
          X_timesp_cov <- cbind(X_psi_yearcov * X_sp, X_psi_yearcov * Y_sp)
          
          numTimeSpaceCov <- numTimeSpaceCov + 2
          
        } else {
          X_timesp_cov <- NULL
        }
        
      }
      
      # standard covariates
      {
        
        column_covariate_psi <- ExtractCovariatesFromText(covariates_psi_text)
        
        nameVariables_psi <- colnames(data)[column_covariate_psi]
        ncov_psi <- length(nameVariables_psi)
        
        if(ncov_psi > 0){
          
          X_psi_cov <- stats::model.matrix(~ ., data = X_psi[,column_covariate_psi,drop = F])[,-1]  
          
        } 
        
      }
      
      X_psi <- cbind(X_psi_year)
      if(usingSpatial) X_psi <- cbind(X_psi, X_timesp_cov)
      if(ncov_psi > 0) X_psi <- cbind(X_psi, X_psi_cov)
      
    }
    
    # covariates for p
    
    {
      data_p <- data
      
      column_covariate_p <- ExtractCovariatesFromText(covariates_p_text)
      
      X_p_cov <- data_p[,column_covariate_p, drop = F]
      
      ncov_p <- ncol(X_p_cov)
      
      if(ncov_p > 0){
        nameVariables_p <- colnames(data_p)[column_covariate_p]  
      } else {
        nameVariables_p <- c()
      }
      
      if(usingYearDetProb){
        
        data_p$Year <- as.factor(data_p$Year)
        X_p_year <- stats::model.matrix( ~ . - 1, data = data_p[,c("Year"),drop = F])
        
        X_y_index_p <- apply(X_p_year, 1, function(x) {which(x != 0)})
        X_y_index_p <- unlist(X_y_index_p)
        
        X_p <- X_p_year
      } else {
        X_p <- matrix(1, nrow = nrow(data_p), ncol = 1)
        X_y_index_p <- rep(1, nrow(data_p))
      }
      
      if(ncov_p > 0){
        X_p <- cbind(X_p, X_p_cov)
      } 
      
      X_p <- as.matrix(X_p)
      
      p_intercepts <- ifelse(usingYearDetProb, Y, 1)
      
    }
    
    # data for GOF
    
    {
      
      k_s_all <- data %>% dplyr::arrange(Site, Year)
      
      trueYearStatistics <- k_s_all %>% dplyr::group_by(Year) %>% 
        dplyr::summarise(Detections = sum(Occ)) %>% dplyr::arrange(Year)
      
      if(usingSpatial){
        
        k_s_all$SitePatch <- findClosestPoint(cbind(k_s_all$X_sp,k_s_all$Y_sp), X_tilde)
        
        trueSpaceStatistics <- k_s_all %>% dplyr::group_by(SitePatch) %>% 
          dplyr::summarise(Detections = sum(Occ)) %>% dplyr::arrange(SitePatch)
        
        k_s_all <- k_s_all %>% dplyr::select(SitePatch, Year, ) %>% dplyr::mutate(Present = 0)
        
      } else {
        
        k_s_all <- k_s_all %>% dplyr::select(Year) %>% dplyr::mutate(Present = 0)
        
        trueSpaceStatistics <- NULL
      }
      
    }
  }
  
  # ASSIGN THE PRIOR
  {
    print("Setting up the priors")
    
    # fixed parameters
    {
      phi_psi <- 2
      phi_p <- 2
      
      a_l_T_0 <- 1
      b_l_T_0 <- 1
      
      a_sigma_T_0 <- 2
      b_sigma_T_0 <- .25
      
      a_l_s_0 <- 1
      b_l_s_0 <- 10
      
      a_sigma_s_0 <- 2
      b_sigma_s_0 <- .25
      
      a_sigma_eps <- 2
      b_sigma_eps <- .25
      
      rho_l_T <- .05
      rho_sigma_T <- .001
      
      rho_l_S <- .01
      rho_sigma_S <- .001
    }
    
    sigma_s <- b_sigma_s_0 / (a_sigma_s_0 - 1)
    sigma_T <- b_sigma_T_0 / (a_sigma_T_0 - 1)
    
    l_T <- a_l_T_0 / b_l_T_0
    l_s <- a_l_s_0 / b_l_s_0
    
    mu_psi <- invLogit(prior_psi)
    mu_p <- invLogit(prior_p)
    
    sigma_eps <- .05
    
    # priors on psi
    
    {
      ncolXpsi <- ncol(X_psi) + X_centers
      colnamesXpsi_centers <- NULL
      if(X_centers > 0){
        colnamesXpsi_centers <- paste0("Center",seq_len(X_centers))  
      }
      colnamesXpsi <- c(colnames(X_psi)[1:Y],colnamesXpsi_centers,colnames(X_psi)[-(1:Y)])
      b_psi <- c(rep(mu_psi,Y), rep(0, ncolXpsi - Y))
      B_psi <- matrix(0, nrow = ncolXpsi, ncol = ncolXpsi)
      
      C <- matrix(0, nrow = Y + X_centers, ncol = Y + X_centers)
      C[1:Y, 1:Y] <- K(1:Y,1:Y, sigma_T^2, l_T) + sigma_psi^2
      if(usingSpatial){
        C[Y + 1:X_centers, Y + 1:X_centers] <- K2(X_tilde,X_tilde, sigma_s^2, l_s)  
      }
      
      B_psi[1:(Y + X_centers), 1:(Y + X_centers)] <- C
      
      B_psi[Y + X_centers + seq_len(numTimeSpaceCov), 
            Y + X_centers + seq_len(numTimeSpaceCov)] <- diag(phi_psi^2, nrow = numTimeSpaceCov)
      
      if(ncov_psi > 0){
        
        C_covs <- diag(phi_psi^2, nrow = ncov_psi)
        
        B_psi[(Y + X_centers) + numTimeSpaceCov + 1:ncov_psi,
              (Y + X_centers) + numTimeSpaceCov + 1:ncov_psi] <- C_covs
        
      } 
      
      inv_B_psi <- solve(B_psi)
    }
    
    # prior on p
    
    {
      b_p <- c(rep(mu_p, p_intercepts), rep(0, ncol(X_p) - p_intercepts))
      B_p <- matrix(0, nrow = ncol(X_p), ncol = ncol(X_p))
      diag(B_p)[1:p_intercepts] <- sigma_p^2
      
      if(ncov_p > 0){
        
        C <- diag(phi_p^2, nrow = ncov_p)
        
        B_p[p_intercepts + 1:ncov_p, p_intercepts + 1:ncov_p] <- C
        
      }  
      
      inv_B_p <- solve(B_p)
      
    }
  }
  
  # RUN VB
  {
    # initialize output
    {
      mu_beta_psi_output <- array(NA, dim = c(nruns , maxit, ncolXpsi),
                               dimnames = list(c(), c(), colnamesXpsi))
      # Sigma_beta_psi_output <- array(NA, dim = c(nruns , maxit, ncol(X_psi), ncol(X_psi)),
                               # dimnames = list(c(), c(), colnames(X_psi), colnames(X_psi)))
      
      mu_beta_p_output <- array(NA, dim = c(nruns , maxit, ncol(X_p)),
                             dimnames = list(c(), c(), colnames(X_p)))
      Sigma_beta_p_output <- array(NA, dim = c(nruns , maxit, ncol(X_p), ncol(X_p)),
                               dimnames = list(c(), c(), colnames(X_p), colnames(X_p)))
      
      indexUniqueSite <- which(!duplicated(k_s$Site))
      namesUniqueSite <- k_s$Site[indexUniqueSite]
      # mu_eps_output <- array(NA, dim = c(nruns , maxit, S),
                             # dimnames = list(c(), c(), namesUniqueSite))
      # sd_eps_output <- array(NA, dim = c(nruns , maxit, S),
                             # dimnames = list(c(), c(), namesUniqueSite))
      
      l_T_output <- matrix(NA, nrow = nruns, ncol = maxit)
      sigma_T_output <- array(NA, dim = c(nruns, maxit))

      l_S_output <- matrix(NA, nrow = nruns, ncol = maxit)
      sigma_S_output <- array(NA, dim = c(nruns, maxit))
      
      ELBO_output <- matrix(NA, nrow = nruns, ncol = maxit)
      # sigma_s_output <- array(NA, dim = c(nchain, niter))
      # sigma_eps_output <- array(NA, dim = c(nchain, niter))
      # 
      # psi_mean_output <- array(NA, dim = c(nchain, niter, Y))
      # indexUniqueSite <- which(!duplicated(k_s$Site))
      # namesUniqueSite <- k_s$Site[indexUniqueSite]
      # 
      # gofYear_output <- array(NA, dim = c(nchain, niter, Y))
      # if(usingSpatial){
      #   gofSpace_output <- array(NA, dim = c(nchain, niter, nrow(X_tilde))) 
      # } else {
      #   gofSpace_output <- NULL
      # }
      # 
      # if(storeRE){
      #   eps_unique_output <- array(NA, dim = c(nchain, niter, S),
      #                              dimnames = list(c(), c(), namesUniqueSite))
      # } else {
      #   eps_unique_output <- matrix(0, nchain , S)
      #   colnames(eps_unique_output) <- namesUniqueSite
      # }
      
      lastIter <- rep(NA, nruns)
      
    }
    
    for (run in 1:nruns) {
      
      # initialize parameters
      {
        
        print("Initializing the parameters")
        
        mu_beta_psi <- b_psi
        Sigma_beta_psi <- B_psi
        # psi <- as.vector(logit(X_psi %*% beta_psi))
        
        # Xbetapsi <- X_psi %*% beta_psi
        
        mu_beta_p <- b_p
        Sigma_beta_p <- B_p
        # p <- as.vector(logit(X_p %*% beta_p))
        
        # Xbetap <- X_p %*% beta_p
        
        mu_eps <- rep(0, nrow(k_s))
        sd_eps <- rep(0.01, nrow(k_s))
        
        gamma_z <- rep(NA, nrow(k_s))
        for (i in 1:nrow(k_s)) {
          if(k_s[i,1] == 1){
            gamma_z[i] <- 1
          } else {
            gamma_z[i] <- rbeta(1, 1, 1)
          }
        }
        
        # presences across each sampling occasion
        # z_all <- z[indexes_occobs]
        
        # hyperparameters
        {
          l_T <- a_l_T_0 / b_l_T_0
          
          a_sigma_T <- a_sigma_T_0 
          b_sigma_T <- b_sigma_T_0 
          
          l_S <- a_l_s_0 / b_l_s_0
          
          a_sigma_S <- a_sigma_s_0
          b_sigma_S <- b_sigma_s_0
          
        }
        
        if(usingSpatial){
          
          length_grid_ls <- 50
          l_s_grid <- seq(0.01, 0.5, length.out = length_grid_ls)
          K_s_grid <- array(NA, dim = c(X_centers, X_centers, length_grid_ls))
          inv_K_s_grid <- array(NA, dim = c(X_centers, X_centers, length_grid_ls))
          diag_K_s_grid <- array(NA, dim = c(X_centers, length_grid_ls))
          inv_chol_K_s_grid <- array(NA, dim = c(X_centers, X_centers, length_grid_ls))
          
          for (j in 1:length_grid_ls) {
            l_s <- l_s_grid[j]
            K_s_grid[,,j] <- K2(X_tilde, X_tilde, 1, l_s) + diag(exp(-10), nrow = nrow(X_tilde))
            diag_K_s_grid[,j] <- FastGP::rcppeigen_get_diag(K_s_grid[,,j])
            inv_K_s_grid[,,j] <- FastGP::rcppeigen_invert_matrix(K_s_grid[,,j])
            inv_chol_K_s_grid[,,j] <- FastGP::rcppeigen_invert_matrix(FastGP::rcppeigen_get_chol(K_s_grid[,,j]))
          }  
          
        }
        
      }
      
      for (iter in seq_len(maxit)) {
        
        if(verbose){
          print(paste0("Run = ",run," - Iteration = ",iter))
          if(iter > 2){
            print(paste0("Y.E. tol. = ", round(mu_beta_psi_year_tol, 5),
                         " / S.E. tol. = ", round(mu_beta_psi_spatial_tol, 5),
                         " / C.E. tol. = ", round(mu_beta_psi_cov_tol, 5),
                         " / Sigma psi tol. = ", round(Sigma_beta_psi_tol, 5),
                         " / mu_p tol. = ", round(mu_beta_p_tol, 5),
                         " / mu_p cov. tol. = ", round(mu_beta_p_cov_tol, 5),
                         " / Sigma p tol. = ", round(Sigma_beta_p_tol, 5),
                         " / Maxtol = ", round(maxtol, 5)))
            # print(paste0("l_S = ",l_S, " / sigma_S = ",mean_sigmasq_S))
          }
        }
        
        # print(paste0("a_sigma_S = ",a_sigma_S, " / b_sigma_S = ",b_sigma_S))
        
        # precomute products ----
        
        tX_psi_mu_beta_psi <- tXmubetapsi(X_psi, X_y_index, X_s_index,
                                          mu_beta_psi, Y, X_centers,
                                          ncov_psi + numTimeSpaceCov)
        
        tX_p_mu_beta_p <- X_p %*% mu_beta_p
        # tX_p_mu_beta_p <- tX_p_mu_beta_p
        
        XSigmapsiX_all <- XSigmapsiX(X_psi, X_y_index, X_s_index, Sigma_beta_psi, Y, X_centers,
                                     ncov_psi + numTimeSpaceCov)
        # XSigmapsiX_all <- X_psi %*% Sigma_beta_psi %*% t(X_psi)
      
        XSigmapX_all <- XSigmapX(X_p, X_y_index_p, Sigma_beta_p, p_intercepts, ncov_p)
        
        # update z ----
        
        # gamma_z <- update_gamma_z_cpp(gamma_z, as.matrix(k_s), mu_eps,
        #                           mu_beta_psi, X_psi,
        #                           mu_beta_p, X_p,
        #                           tX_psi_mu_beta_psi,
        #                           tX_p_mu_beta_p)
        gamma_z <- update_gamma_z_cpp(gamma_z, as.matrix(k_s),
                                             mu_eps, sd_eps,
                                             XSigmapsiX_all, XSigmapX_all,
                                             X_y_index, X_s_index,
                                             Y, X_centers, ncov_psi + numTimeSpaceCov,
                                             X_y_index_p, p_intercepts, ncov_p,
                                             Sigma_beta_psi, Sigma_beta_p,
                                             tX_psi_mu_beta_psi,
                                             tX_p_mu_beta_p)
        
        gamma_z_all <- gamma_z[indexes_occobs]
        
        # sample psi ----
        
        # list_psi <- update_psi_cpp_old(mu_beta_psi, Sigma_beta_psi,
        #                            gamma_z, X_psi, b_psi, inv_B_psi,
        #                           Y, X_centers, ncov_psi,
        #                            numTimeSpaceCov, X_y_index, X_s_index,
        #                            mu_eps, sd_eps, sigma_eps)
        list_psi <- update_psi_cpp(mu_beta_psi, Sigma_beta_psi,
                                   XSigmapsiX_all,
                                   tX_psi_mu_beta_psi,
                                   gamma_z, 
                                   X_psi, b_psi, inv_B_psi,
                                   k_s$Site, sort(unique(k_s$Site)),
                                   Y, X_centers, ncov_psi,
                                   numTimeSpaceCov, X_y_index, X_s_index,
                                   mu_eps, sd_eps, sigma_eps)
        mu_beta_psi <- list_psi$mu_beta_psi
        Sigma_beta_psi <- list_psi$Sigma_beta_psi
        mu_eps <- list_psi$mu_eps
        sd_eps <- list_psi$sd_eps
        
        # sample hyperparameters ---------
        
        # list_l_sigma_T <- update_l_sigma_T(mu_beta_psi[1:Y],
        #                                    l_T,
        #                                    sigma_T,
        #                                    sigma_psi, 1:Y,
        #                                    a_l_T_0, b_l_T_0,
        #                                    a_sigma_T_0, b_sigma_T_0)
        # sigma_T <- list_l_sigma_T$sigma_T
        # l_T <- list_l_sigma_T$l_T
        
        l_T <- update_l_T_new(mu_beta_psi[1:Y],
                              l_T,
                              a_sigma_T, b_sigma_T,
                              sigma_psi, 1:Y,
                              a_l_T_0, b_l_T_0,
                              a_sigma_T_0, b_sigma_T_0)
        
        list_sigma <- update_sigmasq_T_cpp(mu_beta_psi[1:Y],
                                           l_T,
                                           sigma_psi, 1:Y,
                                           a_sigma_T_0, b_sigma_T_0)
        a_sigma_T <- list_sigma$a_sig
        b_sigma_T <- list_sigma$b_sig

              
        # update inv_B_psi
        {
          mean_sigmasq_T <- b_sigma_T / (a_sigma_T - 1)
          mean_l_T <- l_T

          K_l <- K(1:Y, 1:Y, mean_sigmasq_T, mean_l_T) + sigma_psi^2 + diag(exp(-8), nrow = Y)

          inverse_Kl <- FastGP::rcppeigen_invert_matrix(K_l)
          inverse_Kl[lower.tri(K_l)] <- t(inverse_Kl)[lower.tri(K_l)]

          inv_B_psi[1:Y, 1:Y] <- inverse_Kl
        }

        if(usingSpatial){

          # l_S <- update_l_S_new(mu_beta_psi[Y + 1:X_centers],
          #                                    l_S, a_sigma_S, b_sigma_S,
          #                                    0, X_tilde,
          #                                    a_l_s_0, b_l_s_0,
          #                                    a_sigma_s_0, b_sigma_s_0)
          
          list_l_S <- update_l_S_grid(mu_beta_psi[Y + 1:X_centers],
                                l_s_grid, diag_K_s_grid,
                                inv_K_s_grid,
                                a_sigma_S, b_sigma_S,
                                a_l_s_0, b_l_s_0)
          l_S <- list_l_S$l_S
          inv_Kl <- list_l_S$inv_Kl

          list_sigma <- update_sigmasq_S_cpp_fast(mu_beta_psi[Y + 1:X_centers],
                                             l_S, inv_Kl,
                                             a_sigma_s_0, b_sigma_s_0)
          a_sigma_S <- list_sigma$a_sig
          b_sigma_S <- list_sigma$b_sig
          
          # list_sigma <- update_sigmasq_S_cpp(mu_beta_psi[Y + 1:X_centers],
          #                                    l_S,
          #                                    X_tilde, 0,
          #                                    a_sigma_s_0, b_sigma_s_0)
          # a_sigma_S <- list_sigma$a_sig
          # b_sigma_S <- list_sigma$b_sig
          
          # l_S <- update_l_S(mu_beta_psi[Y + 1:X_centers],
          #                        l_S,
          #                        a_sigma_S, b_sigma_S,
          #                        0, X_tilde,
          #                        a_l_s_0, b_l_s_0,
          #                        sims_grad = 30)
          # 
          # 
          # list_sigma_S <- update_sig_S_delta_cpp(mu_beta_psi[Y + 1:X_centers],
          #                              l_S,
          #                              a_sigma_S, b_sigma_S,
          #                              0, X_tilde,
          #                              a_sigma_s_0, b_sigma_s_0,
          #                              sims_grad = 100, rho_sigma_S)
          # a_sigma_S <- list_sigma_S$a_sig
          # b_sigma_S <- list_sigma_S$b_sig
          # 
          # c(l_S, a_sigma_S, b_sigma_S)
          # 
          mean_sigmasq_S <- b_sigma_S / (a_sigma_S - 1)

          # K_l <- K2(X_tilde, X_tilde, mean_sigmasq_S, l_S) +
          #   diag(exp(-8), nrow = nrow(X_tilde))
          # 
          # inverse_Kl <- FastGP::rcppeigen_invert_matrix(K_l)
          # inverse_Kl[lower.tri(K_l)] <- t(inverse_Kl)[lower.tri(K_l)]
          # 
          # inv_B_psi[Y + 1:X_centers, Y + 1:X_centers] <- inverse_Kl
          inv_B_psi[Y + 1:X_centers, Y + 1:X_centers] <- (1 / mean_sigmasq_S) * inv_Kl

        }
        
        # update p ----------------
        
        if(sum(gamma_z) > 0){
          # list_p <- update_p(mu_beta_p, Sigma_beta_p, 
          #                    y, gamma_z_all, X_p, b_p, inv_B_p, Xbeta_p, 
          #                    ncov_p, p_intercepts, X_y_index_p, usingYearDetProb)
          list_p <- update_p_cpp(mu_beta_p, Sigma_beta_p, y, gamma_z_all, tX_p_mu_beta_p, 
                                 XSigmapX_all, X_p, b_p, 
                                 inv_B_p, ncov_p, p_intercepts, X_y_index_p, usingYearDetProb)
          mu_beta_p <- list_p$mu_beta_p
          Sigma_beta_p <- list_p$Sigma_beta_p
        }
        
        # compute ELBO -------
        
        # current_ELBOval <- ELBO(y,
        #                         mu_beta_psi, Sigma_beta_psi,
        #                         mu_beta_p, Sigma_beta_p,
        #                         gamma_z,
        #                         l_T, sigma_T,
        #                         l_S, sigma_S)
        
        # write parameters  -------------------------------------------------------
        
        mu_beta_psi_output[run,iter,] <- mu_beta_psi
        # Sigma_beta_psi_output[run,iter,,] <- Sigma_beta_psi
        
        mu_beta_p_output[run,iter,] <- mu_beta_p
        Sigma_beta_p_output[run,iter,,] <- Sigma_beta_p
        
        # mu_eps_output[run,iter,] <- mu_eps[indexUniqueSite]
        # sd_eps_output[run,iter,] <- sd_eps[indexUniqueSite]
        
        l_T_output[run,iter] <- l_T
        
        # sigma_T_output[run,iter] <- sigma_T
        
        if(usingSpatial){
          
          # sigma_S_output[run,iter] <- sigma_S
          
          l_S_output[run,iter] <- l_S
            
        }
        
        # ELBO_output[run,iter] <- current_ELBOval
        
        # check tolerance ----------
        
        insideTolerance <- F
        
        if(iter > 1){
          
          mu_beta_psi_year_tol <- max(abs((mu_beta_psi_output[run,iter,1:Y] - mu_beta_psi_output[run,iter - 1,1:Y]) / 
                                       mu_beta_psi_output[run,iter - 1,1:Y]))
          
          if(usingSpatial){
            mu_beta_psi_spatial_tol <- max(abs((mu_beta_psi_output[run,iter,Y + 1:X_centers] - 
                                                  mu_beta_psi_output[run,iter - 1, Y + 1:X_centers]) / 
                                                 mu_beta_psi_output[run,iter - 1, Y + 1:X_centers]))  
          } else {
            mu_beta_psi_spatial_tol <- 0
          }
        
          mu_beta_psi_cov_tol <- max(abs((mu_beta_psi_output[run,iter,Y + X_centers + seq_len(numTimeSpaceCov + ncov_psi)] - 
                                                mu_beta_psi_output[run,iter - 1, Y + X_centers + seq_len(numTimeSpaceCov + ncov_psi)]) / 
                                       mu_beta_psi_output[run,iter - 1, Y + X_centers + seq_len(numTimeSpaceCov + ncov_psi)]))
          
          # Sigma_beta_psi_tol <- max(diag(abs((Sigma_beta_psi_output[run,iter,,] - Sigma_beta_psi_output[run,iter - 1,,]) / 
                                          # Sigma_beta_psi_output[run,iter - 1,,])))
          Sigma_beta_psi_tol <- max(diag(abs((Sigma_beta_psi - Sigma_beta_psi_old) / 
                                               Sigma_beta_psi_old)))
          Sigma_beta_psi_old <- Sigma_beta_psi
          
          # p
          
          mu_beta_p_tol <- max(abs((mu_beta_p_output[run,iter,1:p_intercepts] - mu_beta_p_output[run,iter - 1,1:p_intercepts]) / 
                                     mu_beta_p_output[run,iter - 1,1:p_intercepts]))
          if(ncov_p > 0){
            mu_beta_p_cov_tol <- max(abs((mu_beta_p_output[run,iter, p_intercepts + seq_len(ncov_p)] - 
                                            mu_beta_p_output[run,iter - 1, p_intercepts + seq_len(ncov_p)]) / 
                                           mu_beta_p_output[run,iter - 1, p_intercepts + seq_len(ncov_p)]))  
          } else {
            mu_beta_p_cov_tol <- 0
          }
          
          Sigma_beta_p_tol <- max(diag(abs((Sigma_beta_p_output[run,iter,,] - Sigma_beta_p_output[run,iter - 1,,]) / 
                                        Sigma_beta_p_output[run,iter - 1,,])))
          
          # a_l_T_tol <- max(abs((a_l_T_output[run,iter] - a_l_T_output[run,iter - 1]) / 
          #                        a_l_T_output[run,iter - 1]))
          # b_l_T_tol <- max(abs((b_l_T_output[run,iter] - b_l_T_output[run,iter - 1]) / 
          #                        b_l_T_output[run,iter - 1]))
          
          # l_T_tol <- max(abs(((a_l_T_output[run,iter] / b_l_T_output[run,iter]) - (a_l_T_output[run,iter - 1] / b_l_T_output[run,iter - 1])) /
                                 # (a_l_T_output[run,iter - 1] / b_l_T_output[run,iter - 1])))
          
          # sigma_T_tol <- max(abs(((a_sigma_T_output[run,iter] / b_sigma_T_output[run,iter]) - (a_sigma_T_output[run,iter - 1] / b_sigma_T_output[run,iter - 1])) / 
          #                        (a_sigma_T_output[run,iter - 1] / b_sigma_T_output[run,iter - 1])))
          
          maxtol <- max(mu_beta_psi_year_tol, mu_beta_psi_spatial_tol, mu_beta_psi_cov_tol,
                        Sigma_beta_psi_tol, 
                        mu_beta_p_tol, mu_beta_p_cov_tol,
                        Sigma_beta_p_tol)
          
          if(maxtol < tol){
            insideTolerance <- T
          } 
          
        } else {
          
          Sigma_beta_psi_old <- Sigma_beta_psi
          
        }
        
        if(insideTolerance){
          lastIter[run] <- iter
          break
        }
        # mu_beta_psi_output[chain,iter - nburn,] <- mu_beta_psi
        # Sigma_beta_psi_output[chain,iter - nburn,,] <- Sigma_beta_psi
        # mu_beta_p_output[chain,iter - nburn,] <- mu_beta_p
        # Sigma_beta_p_output[chain,iter - nburn,,] <- Sigma_beta_p
        
      }
      
      if(!insideTolerance){
        lastIter[run] <- maxit
      }
      
    }
    
  }
  
  # COMPARE RUNS
  {
    # for (run in 1:nruns) {
      
      # rangeInLine <- T
      # 
      # # mu beta psi
      # mu_beta_psi_all <- mapply(mu_beta_psi_output)
      # 
      # mu_beta_psi_differences <- apply(mu_beta_psi_all, , function(x){
      #   range(x) / max(x)
      # })
      # 
      # if(any(range > .05)){
      #   rangeInLine <- F
      # }
      
    # }
  }
  
  # store summaries of the dataset
  {
    
    dataCharacteristics <- list("Years" = years,
                                "X_tilde" = X_tilde,
                                "originalsites" = originalsites,
                                "gridStep" = gridStep,
                                "indexUniqueSite" = indexUniqueSite,
                                "XY_sp_unique" = XY_sp_unique,
                                "XY_centers_unique" = XY_centers_unique,
                                "usingYearDetProb" = usingYearDetProb,
                                "usingSpatial" = usingSpatial,
                                "numTimeSpaceCov" = numTimeSpaceCov,
                                "X_psi_yearcov_values" = X_psi_yearcov_values,
                                "nameVariables_psi" = nameVariables_psi,
                                "nameVariables_p" = nameVariables_p) 
    
  }
  
  # create model output
  {
    # GOF_output <- list(
    #   "trueYearStatistics" = trueYearStatistics,
    #   "trueSpaceStatistics" = trueSpaceStatistics,
    #   "gofYear_output" = gofYear_output,
    #   "gofSpace_output" = gofSpace_output
    # )
    
    modelOutput <- list(
      "mu_beta_psi" = mu_beta_psi_output[1, lastIter[1],],
      "Sigma_beta_psi" = Sigma_beta_psi,
      "mu_beta_p" = mu_beta_p_output[1, lastIter[1],],
      "Sigma_beta_p" = Sigma_beta_p_output[1, lastIter[1],,],
      "mu_eps" = mu_eps,
      "sd_eps" = sd_eps,
      "l_T" = l_T_output[1, lastIter[1]],
      "l_S" = l_S_output[1, lastIter[1]],
      "sigma_T" = sigma_T_output[1, lastIter[1]],
      "sigma_S" = sigma_S_output[1, lastIter[1]],
      "ELBO" = ELBO_output
    )
  }
  
  # modelResults <- list("dataCharacteristics" = dataCharacteristics,
                       # "modelOutput" = modelOutput)
  list("dataCharacteristics" = dataCharacteristics,
       "modelOutput" = modelOutput)
}