
library("Mcomp")
library("Rlgt")
M3.data <- subset(M3,"yearly")

set.seed(12)

curr_series <- 1

sizeTestSet <- length(M3.data[[curr_series]]$xx)
data.train <- M3.data[[curr_series]]$x

mod <- list()
forecasts <- list()


#fc <- forecast(ets(data.train), h = sizeTestSet)
#plot(fc)


##--------------------------------
##Fit LGT model
#
##plot(forecast(fit.lgt(data.train)))
#
##mod[["lgt"]] <- fit.lgt(data.train)
#
#myMod <- init.lgt()
#
##mod[["lgt"]] <- fit.lgt(data.train, h = sizeTestSet, stanModel=myMod, ncores=4)
#
#options(error=recover)
#
mod[["lgt"]] <- fit.lgt(data.train, model="LGT2", nCores=4, nChains=4,
    control=lgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=2000), 
    verbose=TRUE)
forecasts[["lgt"]] <- forecast(mod[["lgt"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["lgt"]],main=paste(curr_series,'by LGT'))

#--------------------------------
#Fit LGTe model
mod[["lgte"]] <- fit.lgt(data.train, model="LGTe", nCores=4, nChains=4,
		control=lgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=2000), 
		verbose=TRUE)
forecasts[["lgte"]] <- forecast(mod[["lgte"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["lgte"]],main=paste(curr_series,'by LGTe'))

#--------------------------------
#Fit Trend model
mod[["trend"]] <- fit.lgt(data.train, model="Trend", nCores=4, nChains=4,
		control=lgt.control(MAX_NUM_OF_REPEATS=2), 
		verbose=TRUE)
forecasts[["trend"]] <- forecast(mod[["trend"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["trend"]],main=paste(curr_series,'by Trend'))


############################################
M3.data <- subset(M3,"monthly")
curr_series=805
		#sample(length(M3.data),1)
sizeTestSet <- length(M3.data[[curr_series]]$xx)
data.train <- M3.data[[curr_series]]$x


##--------------------------------
##Fit SGT model
mod[["sgt"]] <- fit.lgt(data.train, model="SGT2", nCores=2, nChains=4,
		control=lgt.control(MAX_NUM_OF_REPEATS=2), 
		verbose=TRUE)
forecasts[["sgt"]] <- forecast(mod[["sgt"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["sgt"]],main=paste(curr_series,'by SGT'))

#--------------------------------
#Fit LGTe model
mod[["sgte"]] <- fit.lgt(data.train, model="SGTe", nCores=4, nChains=4,
		control=lgt.control(), 
		verbose=TRUE)
forecasts[["sgte"]] <- forecast(mod[["sgte"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["sgte"]],main=paste(curr_series,'by SGTe'))
