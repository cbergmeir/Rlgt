
library("Mcomp")
library("Rlgt")
library("R.utils")
set.seed(12)




# Test Sdampen
M3.data <- c(subset(M3,"quarterly"), subset(M3,"monthly"))
forecast_sdampen <- list()
error <- 0



for (iter in 1001:1500) {
  data.train <- M3.data[[iter]]$x
  data.test <- M3.data[[iter]]$xx
  sizeTestSet <- length(data.test)
  
  
  rstanmodel <- NULL
  rstanmodel <- evalWithTimeout(fit.dampen(data.train, model="SDampen", nCores=1, nChains=4,
                           control=lgt.control(MAX_NUM_OF_REPEATS=10, NUM_OF_ITER=1000),
                           verbose=FALSE), timeout=1000, onTimeout = "warning")
  if (!is.null(rstanmodel)){  
    forecast_sdampen[[iter]] <- forecast(rstanmodel, h = sizeTestSet, level=c(80, 95, 98))
  }
  else{
    forecast_sdampen[[iter]] <- hw(data.train,seasonal="multiplicative", h = sizeTestSet)
    error <- error + 1
  }
}

#saveRDS(c(forecast_dampen, forecast_sdampen), "result1.rds")

myPath <- "/home/rwinwibowo/Documents/research/output"

saveRDS(forecast_sdampen, file.path(myPath, "sdampen3.rds"))
