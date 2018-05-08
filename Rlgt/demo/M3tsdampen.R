
library("Mcomp")
library("Rlgt")
set.seed(12)


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
  forecast_tsdampen[[iter]] <- forecast(rstanmodel, h = sizeTestSet, level=c(80, 95, 98))$mean
}

myPath <- "/home/ubuntu/Documents/Experiment/"

saveRDS(forecast_sdampen, file.path(myPath, "tsdampen.rds"))