
################################################################################
# DATA LOADING AND FORMATTING
################################################################################

# READ IN RAW DATA + LIBRARIES
library(tidyverse)
landscapeData_gEC_RAW <- read_csv("/Users/benedictcummins-mburu/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/EEB397-project/landscapeDataGammaEC.csv")
landscapeData_gO_RAW <-  read_csv("/Users/benedictcummins-mburu/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/EEB397-project/landscapeDataGammaO.csv")
landscapeData_eEC_RAW <- read_csv("/Users/benedictcummins-mburu/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/EEB397-project/landscapeDataE1EC.csv")
landscapeData_eO_RAW <-  read_csv("/Users/benedictcummins-mburu/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/EEB397-project/landscapeDataE1O.csv")


# CONVERT TO CONSUMER DATA (as opposed to coexistance pattern data)
safe_min <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NaN) 
  return(min(x, na.rm = TRUE))
}
convert_to_species <- function(rawLandscapeData) {
  # initialize return data
  return_data <- data.frame(ParamValue = numeric(),
                            PP.2 =       numeric(),
                            PP.3 =       numeric(),
                            MeanCon.2 =  numeric(),
                            MeanCon.3 =  numeric(),
                            MeanAv.2 =   numeric(),
                            MeanAv.3 =   numeric(),
                            MinCon.2 =   numeric(),
                            MinCon.3 =   numeric(),
                            MinAv.2 =    numeric(),
                            MinAv.3 =    numeric()
  )
  # for each parametert value, get all relevant data
  for (i in 1:length(unique(rawLandscapeData$ParameterValue))) {
    # setup
    this_value <- unique(rawLandscapeData$ParameterValue)[i]
    this_dat <- rawLandscapeData %>% filter(ParameterValue == this_value)
    sum_PP2 <- sum(this_dat$PP[this_dat$DistributionID %in% c(0.50, 1.0)], na.rm = TRUE)
    sum_PP3 <- sum(this_dat$PP[this_dat$DistributionID %in% c(0.75, 1.0)], na.rm = TRUE)
    minthing1 <- this_dat$MinConnectivity[this_dat$DistributionID %in% c(0.50, 1.0)]
    minthing2 <- this_dat$MinConnectivity[this_dat$DistributionID %in% c(0.75, 1.0)]
    minthing3 <- this_dat$MinAvailability[this_dat$DistributionID %in% c(0.50, 1.0)]
    minthing4 <- this_dat$MinAvailability[this_dat$DistributionID %in% c(0.75, 1.0)]
    # create temp df
    temp_df <- data.frame(
      ParamValue = this_value,
      PP.2 =       sum_PP2,
      PP.3 =       sum_PP3,
      MeanCon.2 = sum(this_dat$PP[this_dat$DistributionID %in% c(0.50, 1.0)] * this_dat$MeanConnectivity[this_dat$DistributionID %in% c(0.50, 1.0)], na.rm = TRUE) / sum_PP2,
      MeanCon.3 = sum(this_dat$PP[this_dat$DistributionID %in% c(0.75, 1.0)] * this_dat$MeanConnectivity[this_dat$DistributionID %in% c(0.75, 1.0)], na.rm = TRUE) / sum_PP3,
      MeanAv.2 =  sum(this_dat$PP[this_dat$DistributionID %in% c(0.50, 1.0)] * this_dat$MeanAvailability[this_dat$DistributionID %in% c(0.50, 1.0)], na.rm = TRUE) / sum_PP2,
      MeanAv.3 =  sum(this_dat$PP[this_dat$DistributionID %in% c(0.75, 1.0)] * this_dat$MeanAvailability[this_dat$DistributionID %in% c(0.75, 1.0)], na.rm = TRUE) / sum_PP3,
      MinCon.2 =  safe_min(minthing1),
      MinCon.3 =  safe_min(minthing2),
      MinAv.2 =  safe_min(minthing3),
      MinAv.3 =  safe_min(minthing4)
    )
    return_data <- rbind(return_data, temp_df)
  }
  return(return_data)
}
gEC <- convert_to_species(landscapeData_gEC_RAW)
gO <- convert_to_species(landscapeData_gO_RAW)
eEC <- convert_to_species(landscapeData_eEC_RAW)
eO <- convert_to_species(landscapeData_eO_RAW)

################################################################################
# DATA PLOTTING: Gamma effects in EC
################################################################################

# plot mean connectivity
plot(x = gEC$ParamValue, y = gEC$MeanCon.2,
     ylim = c(0, 1), ylab = "Mean Connectivity", xlab = "Gamma",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = gEC$ParamValue, y = gEC$MeanCon.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")
# plot mean availability
plot(x = gEC$ParamValue, y = gEC$MeanAv.2,
     ylim = c(0, 1), ylab = "Mean Availability", xlab = "Gamma",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = gEC$ParamValue, y = gEC$MeanAv.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")


# plot min connectivity
plot(x = gEC$ParamValue, y = gEC$MinCon.2,
     ylim = c(0, 1), ylab = "Min Connectivity", xlab = "Gamma",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = gEC$ParamValue, y = gEC$MinCon.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")
# plot mean availability
plot(x = gEC$ParamValue, y = gEC$MinAv.2,
     ylim = c(0, 1), ylab = "Min Availability", xlab = "Gamma",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = gEC$ParamValue, y = gEC$MinAv.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")

# plot pp
plot(x = gEC$ParamValue, y = gEC$PP.2,
     ylim = c(0, 1), ylab = "Persistence Extent", xlab = "Gamma",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = gEC$ParamValue, y = gEC$PP.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")


################################################################################
# DATA PLOTTING: Gamma effects in O
################################################################################

# plot mean connectivity
plot(x = gO$ParamValue, y = gO$MeanCon.2,
     ylim = c(0, 1), ylab = "Mean Connectivity", xlab = "Gamma",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = gO$ParamValue, y = gO$MeanCon.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")
# plot mean availability
plot(x = gO$ParamValue, y = gO$MeanAv.2,
     ylim = c(0, 1), ylab = "Mean Availability", xlab = "Gamma",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = gO$ParamValue, y = gO$MeanAv.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")


# plot min connectivity
plot(x = gO$ParamValue, y = gO$MinCon.2,
     ylim = c(0, 1), ylab = "Min Connectivity", xlab = "Gamma",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = gO$ParamValue, y = gO$MinCon.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")
# plot mean availability
plot(x = gO$ParamValue, y = gO$MinAv.2,
     ylim = c(0, 1), ylab = "Min Availability", xlab = "Gamma",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = gO$ParamValue, y = gO$MinAv.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")


# plot pp
plot(x = gO$ParamValue, y = gO$PP.2,
     ylim = c(0, 1), ylab = "Persistence Extent", xlab = "Gamma",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = gO$ParamValue, y = gO$PP.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")





















################################################################################
# DATA PLOTTING: E1 effects in EC
################################################################################

# plot mean connectivity
plot(x = eEC$ParamValue, y = eEC$MeanCon.2,
     ylim = c(0, 1), ylab = "Mean Connectivity", xlab = "Intrinsic Extinction",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = eEC$ParamValue, y = eEC$MeanCon.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")
# plot mean availability
plot(x = eEC$ParamValue, y = eEC$MeanAv.2,
     ylim = c(0, 1), ylab = "Mean Availability", xlab = "Intrinsic Extinction",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = eEC$ParamValue, y = eEC$MeanAv.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")


# plot min connectivity
plot(x = eEC$ParamValue, y = eEC$MinCon.2,
     ylim = c(0, 1), ylab = "Min Connectivity", xlab = "Intrinsic Extinction",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = eEC$ParamValue, y = eEC$MinCon.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")
# plot mean availability
plot(x = eEC$ParamValue, y = eEC$MinAv.2,
     ylim = c(0, 1), ylab = "Min Availability", xlab = "Intrinsic Extinction",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = eEC$ParamValue, y = eEC$MinAv.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")

# plot pp
plot(x = eEC$ParamValue, y = eEC$PP.2,
     ylim = c(0, 1), ylab = "Persistence Extent", xlab = "Intrinsic Extinction",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = eEC$ParamValue, y = eEC$PP.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")


################################################################################
# DATA PLOTTING: E1 effects in O
################################################################################

# plot mean connectivity
plot(x = eO$ParamValue, y = eO$MeanCon.2,
     ylim = c(0, 1), ylab = "Mean Connectivity", xlab = "Intrinsic Extinction",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = eO$ParamValue, y = eO$MeanCon.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")
# plot mean availability
plot(x = eO$ParamValue, y = eO$MeanAv.2,
     ylim = c(0, 1), ylab = "Mean Availability", xlab = "Intrinsic Extinction",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = eO$ParamValue, y = eO$MeanAv.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")


# plot min connectivity
plot(x = eO$ParamValue, y = eO$MinCon.2,
     ylim = c(0, 1), ylab = "Min Connectivity", xlab = "Intrinsic Extinction",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = eO$ParamValue, y = eO$MinCon.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")
# plot min availability
plot(x = eO$ParamValue, y = eO$MinAv.2,
     ylim = c(0, 1), ylab = "Min Availability", xlab = "Intrinsic Extinction",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = eO$ParamValue, y = eO$MinAv.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")


# plot pp
plot(x = eO$ParamValue, y = eO$PP.2,
     ylim = c(0, 1), ylab = "Persistence Extent", xlab = "Intrinsic Extinction",
     type = "o", pch = 22, col = "#21908C", bg = "#21908C")
lines(x = eO$ParamValue, y = eO$PP.3,
      type = "o", pch = 21, col = "#5DC863", bg = "#5DC863")


