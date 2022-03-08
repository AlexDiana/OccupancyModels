library(MASS); library(dplyr); 
library(coda); library(FastGP); library(here)
library(reshape2); library(ggplot2); 
library(Rcpp); library(RcppArmadillo)

sourceCpp(here("MCMC/Functions","codecpp.cpp"))
r_functions <- list.files(here("MCMC/Functions"))
r_functions <- r_functions[grep(".r", r_functions)]
sapply(r_functions, function(x) {
  source(here("MCMC/Functions",x))
})

# true data ----------

realData <- F

# input
{
  
  if(realData)   { # real data 
    
    data <- read.csv(file = "Ringlet_BNM_1970_2014_processed.csv", stringsAsFactors = F)
    # # data <- read.csv(file = "Duke of Burgundy_BNM_1970_2014_processed.csv", stringsAsFactors = F)
    
    {
      # load("/cluster/home/osr/ad625/FastOccupancy/BFdata/finalmaxtemperature12km.rda")
      # data$MaxTemperature <- temperature
      
      # load("/cluster/home/osr/ad625/FastOccupancy/BFdata/final_mintemperature12km.rda")
      # data$MinTemperature <- temperature
      # load("/cluster/home/osr/ad625/FastOccupancy/BFdata/finalrainfall12km.rda")
      # data$LogRainfall <- log(rainfall + 1)
      # data$isNorth <- data$NORTH > 450000
      # load("/cluster/home/osr/ad625/FastOccupancy/BFdata/final_yearlywind5km.rda")
      # data$YearlyWind <- yearlywind
      
      load("C:/Users/Alex/Dropbox/R Folder/PostDoc/VB Occupancy/VB model/Byron Data/maxListLAreaRinglet.rda")
      # load("/cluster/home/osr/ad625/FastOccupancy/BFdata/maxListLAreaDuke.rda")
      data$maxListlArea <- maxListlArea
      data$relativeListLength <- data$listL / data$maxListlArea
    }
    #
    #
    # subset
    {
      # data <- data[data$Year < 2000,]
    }
    meanDataEast <- mean(data$EAST); sdDataEast <- sd(data$EAST)
    meanDataNorth <- mean(data$NORTH); sdDataNorth <- sd(data$NORTH)
    sdBoth <- (sdDataEast + sdDataNorth) / 2
    data$EAST <- (data$EAST - meanDataEast) / sdBoth
    data$NORTH <- (data$NORTH - meanDataNorth) / sdBoth
    # data$listL <- (data$listL - mean(data$listL)) / sd(data$listL)
    data$relativeListLength <- (data$relativeListLength -
                                  mean(data$relativeListLength)) / sd(data$relativeListLength)
    # # data$MaxTemperature <- (data$MaxTemperature - mean(data$MaxTemperature)) / sd(data$MaxTemperature)
    # # data$MinTemperature <- (data$MinTemperature - mean(data$MinTemperature)) / sd(data$MinTemperature)
    # # data$LogRainfall <- (data$LogRainfall - mean(data$LogRainfall)) / sd(data$LogRainfall)
    # # data$isNorth <- (data$isNorth - mean(data$isNorth)) / sd(data$isNorth)
    # # data$YearlyWind <-  (data$YearlyWind - mean(data$YearlyWind)) / sd(data$YearlyWind)
    {
      library(lubridate)
      data_date <- as.Date(data$Date)
      # monthDate <- month(data_date)
      # data$Season <- ifelse(monthDate %in% c(10,11,12,1,2,3), "COLD","WARM")
      # data$SeasonCold <- data$Season == "COLD"
      data_jDate <- as.POSIXlt(data_date, format = "%d%b%y")
      data$JulianDate <- data_jDate$yday
      data$JulianDateSq <- (data_jDate$yday)^2
      data$JulianDateCub <- (data_jDate$yday)^3
      
      
      data$JulianDate <- scale(data$JulianDate)
      data$JulianDateSq <- scale(data$JulianDateSq)
      data$JulianDateCub <- scale(data$JulianDateCub)
                                  }
    
    # data$YearSt <- (data$Year - mean(data$Year)) / sd(data$Year)
    index_year <- 9
    index_site <- 3
    index_occ <- 7
    index_spatial_x <- 4
    index_spatial_y <- 5
    covariates_p_text <- "11,12,13,14"
    covariates_psi_text <- "0"
  } 
  else  # simulated data
  {
    data_file <- here("Data/Simulated","simdata.csv")
    
    data <- read.csv(file = data_file, stringsAsFactors = F)
    index_year <- 1
    index_site <- 2
    index_occ <- 8
    index_spatial_x <- 3
    index_spatial_y <- 4
    covariates_psi_text <- "5-6"#  "7"
    covariates_p_text <- "7"  #"5-6"
  }
  
}

# prior parameters
{
  prior_psi <- .1
  sigma_psi <- 2
  
  prior_p <- .1
  sigma_p <- 2
  
  gridStep <-  .2025#.135#.2#.135
  # gridStep * sdBoth / 1000
  # buildSpatialGrid(data$EAST, data$NORTH, gridStep)
  
  gridStep <- .45#.1015#
  buildSpatialGrid(data$X1, data$X2, gridStep)
  
  # buildSpatialGrid(data$X1, data$X2, .1)
}

# mcmc parameter
{
  nchain <- 1
  nburn <- 2000#15000
  niter <- 2000#20000
}

usingSpatial <- F
usingYearDetProb <- T

st1 <- Sys.time()
modelResults_MCMC <- runModel(data, 
                         index_year, # index of the column with year
                         index_site, # index of the column with site
                         index_occ, # index of the column with captures
                         index_spatial_x, # index of the column with the x coordinate
                         index_spatial_y, # index of the column with the y coordinate
                         covariates_psi_text, # index of numerical covariates for psi
                         covariates_p_text, # index of numerical covariates for p
                         prior_psi, # prior mean of occupancy probability
                         sigma_psi, # just as in the notes
                         prior_p, 
                         sigma_p, 
                         usingYearDetProb,
                         usingSpatial,
                         gridStep, # width of the spatial grid
                         storeRE = T,
                         nchain, # number of chains
                         nburn, # burn-in iterations
                         niter) # number of iterations
et1 <- Sys.time()
et1-st1

# plot results
plotOccupancyIndex(modelResults) + ylim(c(0,1))
plotSpatialSitesEffect(modelResults)
plotOccupancyCovariatesEffect(modelResults)
plotBaseLineCaptureProb(modelResults)
plotDetectionCovariatesEffect(modelResults)

# print results
printOccupancyIndex(modelResults)
printOccupancyCovariatesEffect(modelResults)
printBaseLineCaptureProb(modelResults)
printDetectionCovariatesEffect(modelResults)

# gof
plotGOFYearDetections(modelResults)
plotGOFSpaceDections(modelResults)

# diagnostics
tracePlot_OccupancyIndexRate(modelResults, 9)
tracePlot_OccupancyYearEffect(modelResults, 5)
tracePlot_OccupancySiteEffect(modelResults, 30)
tracePlot_OccupancyCovariate(modelResults, 1)

tracePlot_DetectionIntercept(modelResults, 1)
tracePlot_DetectionCovariate(modelResults, 4)

tracePlot_sigma_eps(modelResults)

tracePlot_sigma_T(modelResults)
tracePlot_l_T(modelResults)

tracePlot_sigma_S(modelResults)
tracePlot_l_S(modelResults)

