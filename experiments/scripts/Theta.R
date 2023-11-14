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
forecasts <- list()

start.time <- Sys.time()
for (i in 1:length(yearly.data)) {
  print(paste("I'm currently working on....", i))
  series <- yearly.data[[i]]
  n <- length(series)
  trainData <- series[1:(n-H)]
  actuals <- series[(n-H+1):n]
  forec <- thetaf(trainData, h = H, level=c(90,98))

  forecasts[[i]] <- forec
  
  # plot(forec, type = "l")
  # xs <- seq(from=length(trainData)+1,to=length(trainData)+ length(actuals))
  # lines(xs,actuals, col=1, type='b',lwd=2)
  
  sMAPE[i] <- mean(abs(forec$mean-actuals)/(forec$mean+actuals))*200
  MASE[i] <- blgt.MASE(forec$mean, actuals, trainData, 1)
}
end.time <- Sys.time()
print(paste("time difference:", end.time-start.time, attr(end.time-start.time, "units")))
print(paste("sMAPE:", mean(sMAPE), ", MASE:", mean(MASE)))
saveRDS(sMAPE, "results/Theta - sMAPE.yearly.rds")
saveRDS(MASE, "results/Theta - MASE.yearly.rds")
# saveRDS(forecasts, "results/Theta - forecasts.yearly.rds")

############################################################################
# monthly series
print("monthly series====")
H <- 18
sMAPE <- rep(0, length(monthly.data))
MASE <- rep(0, length(monthly.data))
forecasts <- list()

start.time <- Sys.time()
for (i in 1:length(monthly.data)) {
  print(paste("I'm currently working on....", i))
  series <- monthly.data[[i]]
  n <- length(series)
  trainData <- series[1:(n-H)]
  actuals <- series[(n-H+1):n]
  forec <- thetaf(ts(trainData, frequency = 12), h = H, level=c(90,98))
  
  forecasts[[i]] <- forec
  # plot(forec, type = "l")
  # xs <- seq(from=length(trainData)+1,to=length(trainData)+ length(actuals))
  # lines(xs,actuals, col=1, type='b',lwd=2)
  
  sMAPE[i] <- mean(abs(forec$mean-actuals)/(forec$mean+actuals))*200
  MASE[i] <- blgt.MASE(forec$mean, actuals, trainData, 1)
}
end.time <- Sys.time()
print(paste("sMAPE:", mean(sMAPE), ", MASE:", mean(MASE)))
print(paste("time difference:", end.time-start.time, attr(end.time-start.time, "units")))
saveRDS(sMAPE, "results/Theta - sMAPE.monthly.rds")
saveRDS(MASE, "results/Theta - MASE.monthly.rds")
# saveRDS(forecasts, "results/Theta - forecasts.monthly.rds")

############################################################################
# quarterly series
print("quarterly series====")
H <- 8
sMAPE <- rep(0, length(quarterly.data))
MASE <- rep(0, length(quarterly.data))
forecasts <- list()

start.time <- Sys.time()
for (i in 1:length(quarterly.data)) {
  print(paste("I'm currently working on....", i))
  series <- quarterly.data[[i]]
  n <- length(series)
  trainData <- series[1:(n-H)]
  actuals <- series[(n-H+1):n]
  forec <- thetaf(ts(trainData, frequency = 4), h = H, level=c(90,98))
  
  forecasts[[i]] <- forec
  # plot(forec, type = "l")
  # xs <- seq(from=length(trainData)+1,to=length(trainData)+ length(actuals))
  # lines(xs,actuals, col=1, type='b',lwd=2)
  
  sMAPE[i] <- mean(abs(forec$mean-actuals)/(forec$mean+actuals))*200
  MASE[i] <- blgt.MASE(forec$mean, actuals, trainData, 1)
}
end.time <- Sys.time()
print(paste("sMAPE:", mean(sMAPE), ", MASE:", mean(MASE)))
print(paste("time difference:", end.time-start.time, attr(end.time-start.time, "units")))
saveRDS(sMAPE, "results/Theta - sMAPE.quarterly.rds")
saveRDS(MASE, "results/Theta - MASE.quarterly.rds")
# saveRDS(forecasts, "results/Theta - forecasts.quarterly.rds")