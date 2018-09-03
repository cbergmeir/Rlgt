library("Mcomp")
library("Rlgt")
M3.data <- subset(M3,"yearly")

options(width=180)
#set.seed(12)
curr_series <- sample(length(M3.data),1)

sizeTestSet <- length(M3.data[[curr_series]]$xx)
data.train <- M3.data[[curr_series]]$x

mod <- list()
forecasts <- list()

#Fit standard LGT model
mod[["lgt"]] <- rlgt(data.train, 
		control=rlgt.control(MAX_NUM_OF_REPEATS=2, NUM_OF_ITER=2000), 
		verbose=TRUE)
forecasts[["lgt"]] <- forecast(mod[["lgt"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["lgt"]],main=paste(curr_series,'by LGT'))


#Fit LGT model with error size proportional to a smoothed innovation (surprise)
mod[["lgte"]] <- rlgt(data.train, error.size.method="innov", 
		control=rlgt.control(MAX_NUM_OF_REPEATS=2, NUM_OF_ITER=3000), 
		verbose=TRUE)
forecasts[["lgte"]] <- forecast(mod[["lgte"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["lgte"]],main=paste(curr_series,'by LGTe'))

#--------------------------------
#Fit Trend model
#y=data.train; model.type="Trend"; xreg = NULL; control=rlgt.control(MAX_NUM_OF_REPEATS=2); nChains=4; nCores=4; addJitter=TRUE; verbose=TRUE
#mod[["trend"]] <- rlgt(data.train, model.type="Trend", nCores=4, nChains=4,
#		control=rlgt.control(MAX_NUM_OF_REPEATS=2), 
#		verbose=TRUE)
#forecasts[["trend"]] <- forecast(mod[["trend"]], h = sizeTestSet, level=c(80, 95, 98))
#plot(forecasts[["trend"]],main=paste(curr_series,'by Trend'))


############################################
M3.data <- subset(M3,"monthly")
curr_series=sample(length(M3.data),1)
sizeTestSet <- length(M3.data[[curr_series]]$xx)
data.train <- M3.data[[curr_series]]$x
actuals <- M3.data[[curr_series]]$xx

#Fit standard SGT model 
mod[["SGT"]] <- rlgt(data.train, 
		control=rlgt.control(), 
		verbose=TRUE)
forecasts[["SGT"]] <- forecast(mod[["SGT"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["SGT"]],main=paste(curr_series,'by SGT'))
lines(actuals, col=1, lwd=1)	


#Fit SGT model using error size proportional to a smoothed innovation (surprise)
mod[["SGTe"]] <- rlgt(data.train, error.size.method="innov", 
		control=rlgt.control(MAX_NUM_OF_REPEATS=2, NUM_OF_ITER=3000), 
		verbose=TRUE)
forecasts[["SGTe"]] <- forecast(mod[["SGTe"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["SGTe"]],main=paste(curr_series,'by SGTe'))
lines(actuals, col=1, lwd=1)	


#Fit SGT model using generalized seasonality
mod[["gSGT"]] <- rlgt(data.train, seasonality.type="generalized", 
		control=rlgt.control(), 
		verbose=TRUE)
forecasts[["gSGT"]] <- forecast(mod[["gSGT"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["gSGT"]],main=paste(curr_series,'by gSGT'))
lines(actuals, col=1, lwd=1)	


#Fit SGT model using error size proportional to a smoothed innovation and generalized seasonality
mod[["gSGTe"]] <- rlgt(data.train, 
		error.size.method="innov", seasonality.type="generalized",
		control=rlgt.control(), 
		verbose=TRUE)
forecasts[["gSGTe"]] <- forecast(mod[["gSGTe"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["gSGTe"]],main=paste(curr_series,'by gSGTe'))
lines(actuals, col=1, lwd=1)	



