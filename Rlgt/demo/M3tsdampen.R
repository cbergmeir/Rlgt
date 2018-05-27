
library("Mcomp")
library("Rlgt")
library("R.utils")
library("fpp2")

# Test TSdampen
M3.data <- c(subset(M3,"quarterly"), subset(M3,"monthly"))
forecast_tsdampen <- list()
error <- 0

for (iter in 1:500) {
  data.train <- M3.data[[iter]]$x
  data.test <- M3.data[[iter]]$xx
  sizeTestSet <- length(data.test)
  
  rstanmodel <- NULL
  rstanmodel <- try(evalWithTimeout(fit.dampen(data.train, model="TSDampen", nCores=1, nChains=4,
                                               control=lgt.control(MAX_NUM_OF_REPEATS=10, NUM_OF_ITER=1000),
                                               verbose=FALSE), timeout=1000, onTimeout = "silent"))
  if (!is.null(rstanmodel) & !inherits(rstanmodel, "try-error")){  
    forecast_tsdampen[[iter]] <- forecast(rstanmodel, h = sizeTestSet, level=c(80, 95, 98))$mean
  }
  else{
    forecast_tsdampen[[iter]] <- hw(data.train,seasonal="multiplicative", h = sizeTestSet)
    error <- error + 1
  }
  
}

myPath <- "/home/ubuntu/Documents/Experiment"

saveRDS(c(forecast_tsdampen, error), file.path(myPath, "tsdampen1.rds"))