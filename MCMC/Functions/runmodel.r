

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
                     usingSpatial = F, 
                     spatialApprox = "SoD", maxPoints = 10,
                     gridStep,
                     storeRE = F, nchain, nburn, niter, 
                     verbose = T, computeGOF = T){
  
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
          
          X_s_index <- X_psi$Site
          
          X_sp <- X_psi$X_sp
          Y_sp <- X_psi$Y_sp
          
          uniqueIndexesSite <- which(!duplicated(X_psi$Site))
          X_sp_unique <- X_sp[uniqueIndexesSite]
          Y_sp_unique <- Y_sp[uniqueIndexesSite]
          X_s <- cbind(X_sp_unique, Y_sp_unique)
          
          # build the grid
          X_tilde <- as.matrix(buildGrid(X_s, gridStep))
          XY_centers <- findClosestPoint(X_s, X_tilde)
          # get rid of unused cells by refinding the grid
          X_tilde <- X_tilde[sort(unique(XY_centers)),]
          X_centers <- nrow(X_tilde)
          
          if(spatialApprox == "SoD"){
            
            # refit
            # XY_centers <- findClosestPoint(cbind(X_sp,Y_sp), X_tilde)
            X_s_centers <- findClosestPoint(X_s, X_tilde)
            
            # X_psi$Center <- as.factor(XY_centers)
            # X_s_index <- XY_centers
            
            # X_psi_s <- stats::model.matrix( ~ . - 1, data = X_psi[,c("Center"),drop = F])
            
          } else if(spatialApprox == "SoR") {
            
            Xs_Xt_dist <- t(apply(X_s, 1, function(x){
              apply(X_tilde, 1, function(y){
                (x[1] - y[1])^2 + (x[2] - y[2])^2
              })
            }))
            
            X_s_centers <- t(apply(Xs_Xt_dist, 1, function(x){
              order(x)[1:maxPoints]
            }))
            
            # which(abs(XY_sp_unique[X_s_m,] -  cbind(X_sp, Y_sp)))
            # View(cbind(XY_sp_unique[X_s_m,] ,  cbind(X_sp, Y_sp)))
            
          }
          
        } else {
          X_s_index <- numeric(0)
          X_centers <- 0
          X_s <- NULL
          X_s_centers <- NULL
          X_tilde <- NULL
          gridStep <- NULL
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
          colnames(X_timesp_cov) <- c("TxX","TxY")
          
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
          
          X_psi_cov <- stats::model.matrix(~ ., data = X_psi[,column_covariate_psi,drop = F])[,-1,drop = F]  
          
        } else {
          
          X_psi_cov <- matrix(0, nrow = nrow(k_s), ncol = 0)
          
        }
        
      }
      
      if(usingSpatial) X_psi_cov <- cbind(X_psi_cov, X_timesp_cov)
      
      ncol_Xpsi <- Y + X_centers + ncov_psi + numTimeSpaceCov
      
      colnames_Xpsi <- paste0("Year",1:Y)
      if(usingSpatial){
        colnames_Xpsi <- c(colnames_Xpsi, paste0("Center",1:X_centers), colnames(X_timesp_cov))
      }
      if(ncov_psi > 0){
        colnames_Xpsi <- c(colnames_Xpsi, colnames(X_psi_cov))
      }
      
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
      
      a_l_T <- 1
      b_l_T <- 1
      
      a_sigma_T <- 2
      b_sigma_T <- .25
      
      a_l_s <- 1
      b_l_s <- 10
      
      a_sigma_s <- 2
      b_sigma_s <- 1
      
      a_sigma_eps <- 2
      b_sigma_eps <- .25
      
      sd_l_T <- .2 
      sd_sigma_T <- .2 
    }
    
    mu_psi <- invLogit(prior_psi)
    mu_p <- invLogit(prior_p)
    
    b_psi <- c(rep(mu_psi,Y), rep(0, X_centers + numTimeSpaceCov + ncov_psi))
    
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
  
  # run MCMC
  {
    # initialize output
    {
      beta_psi_output <- array(NA, dim = c(nchain , niter, ncol_Xpsi),
                               dimnames = list(c(), c(), colnames_Xpsi))
      
      beta_p_output <- array(NA, dim = c(nchain , niter, ncol(X_p)),
                             dimnames = list(c(), c(), colnames(X_p)))
      
      l_T_output <- matrix(NA, nrow = nchain, ncol = niter)
      l_s_output <- matrix(NA, nrow = nchain, ncol = niter)
      sigma_T_output <- array(NA, dim = c(nchain, niter))
      sigma_s_output <- array(NA, dim = c(nchain, niter))
      sigma_eps_output <- array(NA, dim = c(nchain, niter))
      
      psi_mean_output <- array(NA, dim = c(nchain, niter, Y))
      indexUniqueSite <- which(!duplicated(k_s$Site))
      namesUniqueSite <- k_s$Site[indexUniqueSite]
      
      gofYear_output <- array(NA, dim = c(nchain, niter, Y))
      if(usingSpatial){
        gofSpace_output <- array(NA, dim = c(nchain, niter, nrow(X_tilde))) 
      } else {
        gofSpace_output <- NULL
      }
      
      if(storeRE){
        eps_unique_output <- array(NA, dim = c(nchain, niter, S),
                                   dimnames = list(c(), c(), namesUniqueSite))
      } else {
        eps_unique_output <- matrix(0, nchain , S)
        colnames(eps_unique_output) <- namesUniqueSite
      }
      
    }
    
    for (chain in 1:nchain) {
      
      # initialize parameters
      {
        
        print("Initializing the parameters")
        
        beta_psi <- b_psi
        
        beta_p <- b_p
        p <- as.vector(logit(X_p %*% beta_p))
        
        Xbetap <- X_p %*% beta_p
        
        eps_s <- rep(0, nrow(k_s))
        
        # hyperparameters
        {
          sigma_s <- b_sigma_s / (a_sigma_s - 1)
          sigma_T <- b_sigma_T / (a_sigma_T - 1)
          
          l_T <- a_l_T / b_l_T
          l_s <- a_l_s / b_l_s
          
          sigma_eps <- b_sigma_eps / (a_sigma_eps - 1)
          
          # precompute parameters to update l_s
          if(usingSpatial){
            
            length_grid_ls <- 20
            l_s_grid <- seq(0.01, 0.5, length.out = length_grid_ls)
            index_ls <- floor(length_grid_ls / 2)
            l_s <- l_s_grid[index_ls]
            
            # K_s_grid <- array(NA, dim = c(X_centers, X_centers, length_grid_ls))
            # diag_K_s_grid <- array(NA, dim = c(X_centers, length_grid_ls))
            # inv_chol_K_s_grid <- array(NA, dim = c(X_centers, X_centers, length_grid_ls))
            
            inv_K_s_grid <- array(NA, dim = c(X_centers, X_centers, length_grid_ls))
            
            ldet_grid <- rep(NA, length_grid_ls)
            dim1 <- ifelse(spatialApprox == "SoD" | spatialApprox == "SoR", X_centers, S)
            sq_grid <- array(NA, dim = c(dim1, X_centers, length_grid_ls))
            n_p <- X_centers
            
            X_sor_all <- array(NA, dim = c(S, X_centers, length_grid_ls))
            
            for (j in 1:length_grid_ls) {
              
              l_s_current <- l_s_grid[j]
              
              K_l <- K2(X_tilde, X_tilde, 1, l_s_current) + diag(exp(-10), nrow = nrow(X_tilde))
              # K_s_grid[,,j] <- K2(X_tilde, X_tilde, 1, l_s_current) + diag(exp(-10), nrow = nrow(X_tilde))
              # diag_K_s_grid[,j] <- FastGP::rcppeigen_get_diag(K_s_grid[,,j])
              # inv_chol_K_s_grid[,,j] <- FastGP::rcppeigen_invert_matrix(FastGP::rcppeigen_get_chol(K_s_grid[,,j]))
              
              inv_K_s_grid[,,j] <- FastGP::rcppeigen_invert_matrix(K_l)
              
              if(spatialApprox == "SoD" | spatialApprox == "SoR"){
                
                ldet_grid[j] <- sum(log(FastGP::rcppeigen_get_diag(K_l))) * 2
                sq_grid[,,j] <- FastGP::rcppeigen_get_chol(FastGP::rcppeigen_invert_matrix(K_l))
                
                K_staru <- K2(X_s, X_tilde, 1, l_s_current)
                inv_K_uu <- solve(K2(X_tilde, X_tilde, 1, l_s_current))
                X_sor_all[,,j] <- K_staru %*% inv_K_uu
                
                # X_sor_all[,,j] <- apply(X_sor_all[,,j], 1, function(x){
                #   x / max(x)
                # })
                
              } else if(spatialApprox == "SoC"){
                
                if(verbose == T){
                  print(paste0("Precomputing covariance matrix ",j," out of ",length_grid_ls))  
                }
                
                K_staru <- K2(X_s, X_tilde, 1, l_s_current)
                inv_K_uu <- solve(K2(X_tilde, X_tilde, 1, l_s_current) + 
                                    diag(exp(-20), nrow = X_centers))
                X_sor_all[,,j] <- K_staru %*% inv_K_uu
                
                list_params <- ginvsquare_pseudodet(K_staru, K2(X_tilde, X_tilde, 1, l_s_current))
                n_p <- list_params$nval
                ldet_grid[j] <- list_params$lpdet
                sq_grid[,,j] <- list_params$sq
                
                # K_biggest <- apply(X_sor_all[,,i], 1, function(x){
                #   -sort(-x)[1:maxPoints]
                # })
                
                # which_K_biggest <- apply(X_sor_all[,,i], 1, function(x){
                #   order(-x)[1:maxPoints]
                # })
                # 
                # X_s_indexes_all[,,i] <- t(which_K_biggest)
                # for (j in 1:nrow(X_psi)) {
                #   K_staruKu_all[j,-X_s_index_all[j,,i],i] <- 0
                #   K_staruKu_all[j,,i] <- K_staruKu_all[j,,i] / sum(K_staruKu_all[j,,i])
                # }
                
              }
              
            }  
            
            listparams_spatial <- list("inv_K_s_grid" = inv_K_s_grid,
                                       "ldet_grid" = ldet_grid,
                                       "sq_grid" = sq_grid,
                                       "n_p" = n_p,
                                       "X_sor_all" = X_sor_all)
            
          } else {
            l_s_grid <- NULL
            listparams_spatial <- NULL
          }
          
        }
        
        # priors on psi
        
        {
          
          B_psi <- matrix(0, nrow = Y + X_centers + ncov_psi + numTimeSpaceCov, 
                          ncol = Y + X_centers + ncov_psi + numTimeSpaceCov)
          
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
        
        # a_s <- X_s %*% beta_psi
        
        Xbetapsi <- rep(b_psi[1], nrow(k_s))
        psi <- as.vector(logit(Xbetapsi))
        
        z <- rep(NA, nrow(k_s))
        for (i in 1:nrow(k_s)) {
          if(k_s[i,1] == 1){
            z[i] <- 1
          } else {
            z[i] <- stats::rbinom(1, 1, psi[i])
          }
        }
        
        # presences across each sampling occasion
        z_all <- z[indexes_occobs]
        
        
      }
      
      for (iter in seq_len(nburn + niter)) {
        
        if(verbose){
          if(iter <= nburn){
            print(paste0("Chain = ",chain," - Burn-in Iteration = ",iter))  
          } else {
            print(paste0("Chain = ",chain," - Iteration = ",iter - nburn))
          }  
        }
        # print(paste0("sigma = ",sigma_s," - l_s = ",l_s))
        # # print(paste0("l_s = ",l_s," / sigma_s = ",sigma_s,
        # #              " / max(beta_psi) = ",max(abs(beta_psi[Y  + 1:X_centers]))))
        # 
        # sample z ----
        
        z <- sample_z_cpp(psi, p, as.matrix(k_s[,3:4]))
        z_all <- z[indexes_occobs]
        
        # sample psi ----
        
        if (usingSpatial) {
          
          if(spatialApprox == "SoD"){
            
            X_sc_index <- X_s_centers[X_s_index]
            X_s_sor <- NULL
            
          } else if(spatialApprox == "SoR"){
            
            X_sc_index <-  X_s_centers[X_s_index,]
            X_s_sor <- X_sor_all[X_s_index,,index_ls]
            
          }
          
        } else {
          
          X_sc_index <- X_s_index
          X_s_sor <- NULL
          
        }
        
        list_psi <- update_psi(beta_psi, X_psi_cov, Xbetapsi,
                               b_psi, inv_B_psi, 
                               z, k_s, sites, Y, X_centers, ncov_psi, 
                               X_y_index, X_sc_index, X_s_sor,
                               numTimeSpaceCov,
                               usingSpatial, spatialApprox, eps_s, sigma_eps)
        psi <- list_psi$psi
        beta_psi <- list_psi$beta
        Omega <- list_psi$Omega
        eps_s <- list_psi$eps_s
        XbetaY <- list_psi$XbetaY
        Xbetas <- list_psi$Xbetas
        Xbeta_cov <- list_psi$Xbeta_cov
        Xbetapsi <- list_psi$Xbeta
        
        # sample hyperparameters ---------
        
        list_hyperparameters <- update_hyperparameters(l_T, a_l_T, b_l_T, sd_l_T, 
                                                       sigma_T, a_sigma_T, b_sigma_T, sd_sigma_T,
                                                       Y,
                                                       beta_psi, inv_B_psi, 
                                                       b_psi, sigma_psi,
                                                       l_s_grid, 
                                                       a_l_s, b_l_s, 
                                                       sigma_s, a_sigma_s, b_sigma_s, X_tilde,
                                                       a_sigma_eps, b_sigma_eps,
                                                       usingSpatial, 
                                                       XbetaY, Xbetas, Xbeta_cov,
                                                       eps_s, k_s, 
                                                       z, Omega, X_y_index,
                                                       listparams_spatial,
                                                       spatialApprox)
        l_T <- list_hyperparameters$l_T
        sigma_T <- list_hyperparameters$sigma_T
        sigma_s <- list_hyperparameters$sigma_s
        l_s <- list_hyperparameters$l_s
        index_l_s <- list_hyperparameters$index_l_s
        inv_B_psi <- list_hyperparameters$inv_B_psi
        sigma_eps <- list_hyperparameters$sigma_eps
        
        # sample p ----------------
        
        if(sum(z) > 0){
          list_p <- update_p(beta_p, Occs, z_all, X_p, b_p, inv_B_p, Xbetap,
                             ncov_p, p_intercepts, X_y_index_p, usingYearDetProb)
          p <- list_p$p
          beta_p <- list_p$beta
          Xbetap <- list_p$Xbeta
        }
        
        # write parameters  -------------------------------------------------------
        
        if(iter > nburn){
          
          beta_psi_output[chain,iter - nburn,] <- beta_psi
          
          eps_unique <- eps_s[indexUniqueSite]
          
          if(usingSpatial){
            a_s_unique <- Xbetas[indexUniqueSite] + 
              eps_unique
            
            psi_mean_output[chain,iter - nburn,] <- computeYearEffect(Y, a_s_unique, beta_psi)  
          } else {
            psi_mean_output[chain,iter - nburn,] <- logit(beta_psi[1:Y]) + mean(eps_unique)
          }
          
          if(storeRE){
            if(usingSpatial){
              eps_unique_output[chain, iter - nburn,] <- Xbetas[indexUniqueSite] + eps_s[indexUniqueSite]  
            } else {
              eps_unique_output[chain, iter - nburn,] <- eps_s[indexUniqueSite]  
            }
            
          } else {
            eps_unique_output[chain,] <- eps_unique_output[chain,] + 
              (1 / niter)  * eps_s[indexUniqueSite] 
          }
          
          sigma_eps_output[chain, iter - nburn] <- sigma_eps
          
          l_T_output[chain, iter - nburn] <- l_T
          sigma_T_output[chain, iter - nburn] <- sigma_T
          
          l_s_output[chain, iter - nburn] <- l_s
          sigma_s_output[chain, iter - nburn] <- sigma_s
          
          beta_p_output[chain,iter - nburn,] <- beta_p
          
          # GOF
          if (computeGOF) {
            k_s_all$Present <- simulateDetections(p, z_all)
            
            detectionsByYear <- k_s_all %>% dplyr::group_by(Year) %>% 
              dplyr::summarise(Detections = sum(Present)) 
            
            gofYear_output[chain,iter - nburn,] <- detectionsByYear$Detections
            
            if(usingSpatial){
              
              detectionsByPatch <- k_s_all %>% dplyr::group_by(SitePatch) %>% 
                dplyr::summarise(Detections = sum(Present)) 
              
              gofSpace_output[chain,iter - nburn,] <- detectionsByPatch$Detections
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  # store summaries of the dataset
  {
    
    dataCharacteristics <- list("Years" = years,
                                "X_tilde" = X_tilde,
                                "originalsites" = originalsites,
                                "gridStep" = gridStep,
                                "X_s" = X_s,
                                "usingYearDetProb" = usingYearDetProb,
                                "usingSpatial" = usingSpatial,
                                "numTimeSpaceCov" = numTimeSpaceCov,
                                "X_psi_yearcov_values" = X_psi_yearcov_values,
                                "nameVariables_psi" = nameVariables_psi,
                                "nameVariables_p" = nameVariables_p) 
  }
  
  # create model output
  {
    GOF_output <- list(
      "trueYearStatistics" = trueYearStatistics,
      "trueSpaceStatistics" = trueSpaceStatistics,
      "gofYear_output" = gofYear_output,
      "gofSpace_output" = gofSpace_output
    )
    
    modelOutput <- list(
      "beta_psi_output" = beta_psi_output,
      "eps_unique_output" = eps_unique_output,
      "beta_p_output" = beta_p_output,
      "sigma_T_output" = sigma_T_output,
      "l_T_output" = l_T_output,
      "sigma_s_output" = sigma_s_output,
      "l_s_output" = l_s_output,
      "sigma_eps_output" = sigma_eps_output,
      "psi_mean_output" = psi_mean_output,
      "GOF_output" = GOF_output
    )
  }
  
  list("dataCharacteristics" = dataCharacteristics,
       "modelOutput" = modelOutput)
}