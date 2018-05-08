
library("Rlgt")
set.seed(12)

# Use lynx WWWusage for an example
curr_series <- "WWWusage"
data.train <- WWWusage
sizeTestSet <- length(data.train )


mod <- list()
forecasts <- list()

#--------------------------------
#Fit Dampen model
mod[["Dampen"]] <- fit.dampen(data.train, model="Dampen", nCores=4, nChains=4,
  control=lgt.control(MAX_NUM_OF_REPEATS=10, NUM_OF_ITER=2000), 
  verbose=TRUE)
# print the model details
print(mod[["Dampen"]])

# print the interval for all vars
posterior_interval(mod[["Dampen"]])

forecasts[["Dampen"]] <- forecast(mod[["Dampen"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["Dampen"]],main=paste(curr_series,'by Dampen'))

############### Fit TDampen model ############
mod[["TDampen"]] <- fit.dampen(data.train, model="TDampen", nCores=4, nChains=4,
  control=lgt.control(MAX_NUM_OF_REPEATS=10, NUM_OF_ITER=2000,  MIN_NU=15,
    MAX_NU=100), 
  verbose=TRUE)
# print the model details
print(mod[["TDampen"]])

# print the interval for all vars
posterior_interval(mod[["TDampen"]])

forecasts[["TDampen"]] <- forecast(mod[["TDampen"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["TDampen"]],main=paste(curr_series,'by TDampen'))

# Use AirPassenger data as an example for a seasonal dataset
seasonal_data <- AirPassengers
curr_series <- "AirPassengers"
data.train <- seasonal_data 
sizeTestSet <- frequency(AirPassengers)


#--------------------------------
#Fit sDampen model
mod[["SDampen"]] <- fit.dampen(data.train, model="SDampen", nCores=4, nChains=4,
  control=lgt.control(MAX_NUM_OF_REPEATS=10, NUM_OF_ITER=2000), 
  verbose=TRUE)
# print the model details
print(mod[["SDampen"]])

# print the interval for all vars
posterior_interval(mod[["SDampen"]])

#####################



# Use AirPassanger data as an example for a seasonal dataset
seasonal_data <- AirPassengers
curr_series <- "AirPassengers"
data.train <- seasonal_data 
sizeTestSet <- frequency(AirPassengers)


#--------------------------------
#Fit sDampen model
mod[["SDampen"]] <- fit.dampen(data.train, model="SDampen", nCores=4, nChains=4,
  control=lgt.control(MAX_NUM_OF_REPEATS=10, NUM_OF_ITER=2000), 
  verbose=TRUE)
# print the model details
print(mod[["SDampen"]])

# print the interval for all vars
posterior_interval(mod[["SDampen"]])

forecasts[["SDampen"]] <- forecast(mod[["SDampen"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["SDampen"]],main=paste(curr_series,'by SDampen'))

#--------------------------------
#Fit TsDampen model
mod[["TSDampen"]] <- fit.dampen(data.train, model="TSDampen", nCores=4, nChains=4,
  control=lgt.control(MAX_NUM_OF_REPEATS=10, NUM_OF_ITER=2000), 
  verbose=TRUE)
# print the model details
print(mod[["TSDampen"]])

# print the interval for all vars
posterior_interval(mod[["TSDampen"]])

forecasts[["TSDampen"]] <- forecast(mod[["TSDampen"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["TSDampen"]],main=paste(curr_series,'by TSDampen'))



