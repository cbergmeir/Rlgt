
library("Mcomp")
library("Rlgt")
set.seed(12)


# Test tDampen
M3.data <- c(subset(M3,"yearly"), subset(M3,"other"))
forecast_tdampen <- list()


iter<-1

for (iter in 1:length(M3.data)) {
  data.train <- M3.data[[iter]]$x
  data.test <- M3.data[[iter]]$xx
  rstanmodel <- fit.dampen(data.train, model="TDampen", nCores=4, nChains=4,
                           control=lgt.control(MAX_NUM_OF_REPEATS=10, NUM_OF_ITER=2000), 
                           verbose=TRUE)
  sizeTestSet <- length(data.test)
  forecast_tdampen[[iter]] <- forecast(rstanmodel, h = sizeTestSet, level=c(80, 95, 98))
}

# Test TSdampen
M3.data <- c(subset(M3,"quarterly"), subset(M3,"monthly"))
forecast_tsdampen <- list()


iter<-1

for (iter in 1:length(M3.data)) {
  data.train <- M3.data[[iter]]$x
  data.test <- M3.data[[iter]]$xx
  rstanmodel <- fit.dampen(data.train, model="TSDampen", nCores=4, nChains=4,
                           control=lgt.control(MAX_NUM_OF_REPEATS=10, NUM_OF_ITER=2000), 
                           verbose=TRUE)
  sizeTestSet <- length(data.test)
  forecast_tsdampen[[iter]] <- forecast(rstanmodel, h = sizeTestSet, level=c(80, 95, 98))
}

saveRDS(c(forecast_tdampen, forecast_tsdampen), result2.rds)