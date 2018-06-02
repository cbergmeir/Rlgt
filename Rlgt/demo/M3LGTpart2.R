
library("Mcomp")
library("Rlgt")
set.seed(12)


# Test lgt
M3.data <- c(subset(M3,"yearly"), subset(M3,"other"))
forecast_lgt <- list()


iter<-1

for (iter in 401:length(M3.data)) {
  data.train <- M3.data[[iter]]$x
  data.test <- M3.data[[iter]]$xx
  rstanmodel <- fit.lgt(data.train, model="LGT", nCores=1, nChains=4,
                        control=lgt.control(MAX_NUM_OF_REPEATS=10, NUM_OF_ITER=1000), 
                        verbose=FALSE)
  sizeTestSet <- length(data.test)
  forecast_lgt[[iter]] <- forecast(rstanmodel, h = sizeTestSet, level=c(80, 95, 98))$mean
}



myPath <- "/home/ubuntu/Documents/Experiment"

saveRDS(forecast_lgt, file.path(myPath, "lgt4.rds"))