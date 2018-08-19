
library("Mcomp")
library("Rlgt")
M3.data <- subset(M3,"yearly")

set.seed(12)

curr_series <- 1

sizeTestSet <- length(M3.data[[curr_series]]$xx)
data.train <- M3.data[[curr_series]]$x

mod <- list()
forecasts <- list()


#--------------------------------
#Fit LGTe model
mod[["lgte"]] <- rlgt(data.train, model.type="LGTe", nCores=4, nChains=4,
		control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=2000), 
		verbose=TRUE)
forecasts[["lgte"]] <- forecast(mod[["lgte"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["lgte"]],main=paste(curr_series,'by LGTe'))

#--------------------------------
#Fit Trend model
mod[["trend"]] <- rlgt(data.train, model.type="Trend", nCores=4, nChains=4,
		control=rlgt.control(MAX_NUM_OF_REPEATS=2), 
		verbose=TRUE)
forecasts[["trend"]] <- forecast(mod[["trend"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["trend"]],main=paste(curr_series,'by Trend'))


############################################
M3.data <- subset(M3,"monthly")
curr_series=805
		#sample(length(M3.data),1)
sizeTestSet <- length(M3.data[[curr_series]]$xx)
data.train <- M3.data[[curr_series]]$x


#--------------------------------
#Fit SGTe model
mod[["sgte"]] <- rlgt(data.train, model.type="SGTe", nCores=4, nChains=4,
		control=rlgt.control(), 
		verbose=TRUE)
forecasts[["sgte"]] <- forecast(mod[["sgte"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["sgte"]],main=paste(curr_series,'by SGTe'))
