
library("Mcomp")
library("Rlgt")
set.seed(12)


# Test lgt
M3.data <- c(subset(M3,"yearly"), subset(M3,"other"))
forecast_lgt <- list()


iter<-1

for (iter in 172:400) {
  data.train <- M3.data[[iter]]$x
  data.test <- M3.data[[iter]]$xx
  rstanmodel <- fit.lgt(data.train, model="LGT", nCores=2, nChains=4,
                          control=lgt.control(MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=1000), 
                          verbose=FALSE)
  sizeTestSet <- length(data.test)
  forecast_lgt[[iter]] <- forecast(rstanmodel, h = sizeTestSet, level=c(80, 95, 98))$mean
}



myPath <- "/home/rwinwibowo/Documents/Experiments"

saveRDS(forecast_lgt, file.path(myPath, "lgt3a.rds"))