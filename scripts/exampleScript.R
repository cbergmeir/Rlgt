
library("Mcomp")
library("Rlgt2")
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
mod[["lgt"]] <- fit.lgt(data.train, model="LGT", ncores=4, 
    control=lgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=500), 
    verbose=TRUE)
forecasts[["lgt"]] <- forecast(mod[["lgt"]], h = sizeTestSet)
plot(forecasts[["lgt"]])

#--------------------------------
#Fit Trend model

#mod[["trend"]] <- fit.lgt(data.train, model="Trend", ncores=4)
#forecasts[["trend"]] <- forecast(mod[["trend"]], h = sizeTestSet)
#
#forecasts

#-----------------------------------------
#-----------------------------------------
y <- data.train
h <- 6
 stanModel <- init.lgt()
 control=lgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=500, NUM_OF_TRIALS=500)
 ncores=4
 addJitter=TRUE
 verbose=TRUE
