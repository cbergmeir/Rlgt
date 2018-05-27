
library("Mcomp")
library("Rlgt")
set.seed(12)




# Test Sdampen
M3.data <- c(subset(M3,"quarterly"), subset(M3,"monthly"))
forecast_sdampen <- list()


 iter<-1

 for (iter in 1:length(M3.data)) {
   data.train <- M3.data[[iter]]$x
   data.test <- M3.data[[iter]]$xx
   rstanmodel <- fit.dampen(data.train, model="SDampen", nCores=1, nChains=4,
     control=lgt.control(MAX_NUM_OF_REPEATS=10, NUM_OF_ITER=1000),
     verbose=FALSE)
   sizeTestSet <- length(data.test)
   forecast_sdampen[[iter]] <- forecast(rstanmodel, h = sizeTestSet, level=c(80, 95, 98))
}

#saveRDS(c(forecast_dampen, forecast_sdampen), "result1.rds")

myPath <- "/home/ubuntu/Documents/Experiment"

saveRDS(forecast_sdampen, file.path(myPath, "sdampen.rds"))
