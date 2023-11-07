library(tsdl)
library(forecast)
library(Rlgt)
library(doParallel)
# tsdl
dataset <- tsdl

yearly.data <- list()
quarterly.data <- list()
monthly.data <- list()
for (i in 1:length(dataset)) {
  if (!inherits(dataset[[i]],'mts')) {
    series <- as.numeric(dataset[[i]])
    
    # positive forecasts; remove series has NAs
    if (length(which(series <= 0)) == 0 & length(which(is.na(series))) == 0) {
      # plot(dataset[[i]])
      freq <- frequency(dataset[[i]])
      
      if (freq == 1) {
        if (length(series) > 6)
          yearly.data[[length(yearly.data)+1]] <- series
      } else if (freq == 4) {
        if (length(series) > 8)
          quarterly.data[[length(quarterly.data)+1]] <- series
      } else if (freq == 12) {
        if (length(series) > 18)
          monthly.data[[length(monthly.data)+1]] <- series
      }
    }
  }
}

blgt.MASE <- function(yp, yt, train, m) {
  mae <- mean( abs(yp - yt) )
  
  n <- length(train)
  sNaive <- mean( abs( train[1:(n-m)] - train[(m+1):n] ))
  
  mae / sNaive
}

############################################################################
# yearly series
# non-seasonal series
H <- 6
sMAPE <- rep(0, length(yearly.data))
MASE <- rep(0, length(yearly.data))

start.time <- Sys.time()
for (i in 1:length(yearly.data)) {
  series <- yearly.data[[i]]
  n <- length(series)
  trainData <- series[1:(n-H)]
  actuals <- series[(n-H+1):n]
  rstanmodel <- rlgt(trainData,
                     control=rlgt.control(NUM_OF_ITER=4000),
                     verbose=FALSE)

  forec <- forecast(rstanmodel, h = H, level=c(90,98))

  # plot(forec, type = "l")
  # xs <- seq(from=length(trainData)+1,to=length(trainData)+ length(actuals))
  # lines(xs,actuals, col=1, type='b',lwd=2)

  sMAPE[i] <- mean(abs(forec$mean-actuals)/(forec$mean+actuals))*200
  MASE[i] <- blgt.MASE(forec$mean, actuals, trainData, 1)
}
end.time <- Sys.time()
print(paste("time difference:", end.time-start.time))
print(paste("sMAPE:", mean(sMAPE), ", MASE:", mean(MASE)))

