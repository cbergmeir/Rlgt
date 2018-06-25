
# Use lynx WWWusage for an example
curr_series <- "WWWusage"
data.train <- WWWusage
sizeTestSet <- length(data.train )


mod <- list()
forecasts <- list()

#--------------------------------
#Fit SGT model
mod <- holt(data.train, damped=TRUE, h=sizeTestSet)




plot(mod,main=paste(curr_series,'by fpp2'))




# Use AirPassanger data as an example for a seasonal dataset
seasonal_data <- AirPassengers
curr_series <- "AirPassengers"
data.train <- seasonal_data 
sizeTestSet <- frequency(AirPassengers)


#--------------------------------
#Fit SGT model
mod <- hw(data.train,seasonal="multiplicative", h = sizeTestSet)




plot(mod,main=paste(curr_series,'by fpp2'))

