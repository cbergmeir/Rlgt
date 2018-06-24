
library("Mcomp")
library("Rlgt")
set.seed(12)


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
