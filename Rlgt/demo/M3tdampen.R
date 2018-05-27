library("Mcomp")
library("Rlgt")
library("R.utils")
library("fpp2")

set.seed(12)


# Test tDampen
M3.data <- c(subset(M3,"yearly"), subset(M3,"other"))
forecast_tdampen <- list()
error <- 0


#length(M3.data)
for (iter in 1:length(M3.data)) {
  data.train <- M3.data[[iter]]$x
  data.test <- M3.data[[iter]]$xx
  sizeTestSet <- length(data.test)
  
  rstanmodel <- NULL
  rstanmodel <- try(evalWithTimeout(fit.dampen(data.train, model="TDampen", nCores=1, nChains=4,
                                               control=lgt.control(MAX_NUM_OF_REPEATS=10, NUM_OF_ITER=1000),
                                               verbose=FALSE), timeout=1000, onTimeout = "silent"))
  if (!is.null(rstanmodel) & !inherits(rstanmodel, "try-error")){  
    forecast_tdampen[[iter]] <- forecast(rstanmodel, h = sizeTestSet, level=c(80, 95, 98))$mean
  }
  else{
    forecast_tdampen[[iter]] <- holt(data.train, damped=TRUE, phi = 0.9, h=sizeTestSet)
    error <- error + 1
  }

}


myPath <- "/home/ubuntu/Documents/Experiment/"

saveRDS(c(forecast_tdampen, error), file.path(myPath, "tdampen.rds"))