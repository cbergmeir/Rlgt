library("Rlgt")
#set.seed(12)
options(width=180)

# Use lynx WWWusage for an example
curr_series <- "WWWusage"
data.train <- WWWusage
sizeTestSet <- length(data.train )


mod <- list()
forecasts <- list()
#--------------------------------
#Fit LGT model
mod[["LGT"]] <- rlgt(data.train, model.type ="LGT", nCores=4, nChains=4,
  control=lgt.control(MAX_NUM_OF_REPEATS=10, NUM_OF_ITER=2000), 
  verbose=TRUE)
# print the model details
print(mod[["LGT"]])

# print the interval for all vars
posterior_interval(mod[["LGT"]])

forecasts[["LGT"]] <- forecast(mod[["LGT"]], h = sizeTestSet/2, 
                               level=c(80, 95, 98))
plot(forecasts[["LGT"]], main=paste(curr_series,'by LGT'))

# Use AirPassanger data as an example for a seasonal dataset
seasonal_data <- AirPassengers
curr_series <- "AirPassengers"
data.train <- seasonal_data 
sizeTestSet <- frequency(AirPassengers)
#--------------------------------
#Fit SGT model
mod[["SGT"]] <- rlgt(data.train, model.type="SGT", nCores=4, nChains=4,
  control=lgt.control(MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=1000), 
  verbose=TRUE)
# print the model details
print(mod[["SGT"]])

# print the interval for all vars
posterior_interval(mod[["SGT"]])

forecasts[["SGT"]] <- forecast(mod[["SGT"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["SGT"]],main=paste(curr_series,'by SGT'))



