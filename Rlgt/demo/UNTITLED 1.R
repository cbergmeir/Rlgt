
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


for (iter in 401:length(M3.data)) {
  data.train <- M3.data[[iter]]$x
  data.test <- M3.data[[iter]]$xx
  rstanmodel <- fit.lgt(data.train, model="LGT", nCores=1, nChains=4,
                        control=lgt.control(MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=1000), 
                        verbose=FALSE)
  sizeTestSet <- length(data.test)
  forecast_lgt[[iter]] <- forecast(rstanmodel, h = sizeTestSet, level=c(80, 95, 98))$mean
}

for (k in 646:2829){
  data.test[[k]] <- M3[[k]]$xx
  
  M3Forecast[["DAMPEN"]][k]
}
