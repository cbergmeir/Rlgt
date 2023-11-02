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

############################################################################
# yearly series
# non-seasonal series
H <- 6
train.data = list()
future.data = list()
for (i in 1:length(yearly.data))
{
  n <- length(yearly.data[[i]])
  train.data[[i]] = as.numeric(yearly.data[[i]][1:(n-H)])
  future.data[[i]] = as.numeric(yearly.data[[i]][(n-H+1):n])
}
w.series = 1:length(yearly.data)

# homo error
print("yearly tsdl homo error.....")
s = system.time({
  rv=Rlgt::blgt.multi.forecast(train.data[w.series], future.data[w.series], n.samples=1e4, homoscedastic = T)
})
print(s)                         # overall timing info
s[[3]] / length(w.series) # per series time
print("performance in terms of mean sMAPE")
print(mean(rv$sMAPE))
print("performance in terms of mean MASE")
print(mean(rv$MASE))
print("coverage of prediction intervals -- should be close to 95%")
print(mean(rv$InCI)/H)

saveRDS(rv, "tsdl.yearly.homo.rds")

# hetero error
print("yearly tsdl hetero error.....")
s = system.time({
  rv=Rlgt::blgt.multi.forecast(train.data[w.series], future.data[w.series], n.samples=1e4, homoscedastic = F)
})
print(s)                         # overall timing info
s[[3]] / length(w.series) # per series time
print("performance in terms of mean sMAPE")
print(mean(rv$sMAPE))
print("performance in terms of mean MASE")
print(mean(rv$MASE))
print("coverage of prediction intervals -- should be close to 95%")
print(mean(rv$InCI)/H)

saveRDS(rv, "tsdl.yearly.hetero.rds")

############################################################################
# monthly series
H <- 18
train.data = list()
future.data = list()
for (i in 1:length(monthly.data))
{
  n <- length(monthly.data[[i]])
  train.data[[i]] = as.numeric(monthly.data[[i]][1:(n-H)])
  future.data[[i]] = as.numeric(monthly.data[[i]][(n-H+1):n])
}
w.series = 1:length(monthly.data)

# homo error
print("monthly tsdl homo error.....")
s = system.time({
  rv=Rlgt::blgt.multi.forecast(train.data[w.series], future.data[w.series], n.samples=1e4, m = 12, homoscedastic = T)
})
print(s)                         # overall timing info
s[[3]] / length(w.series) # per series time
print("performance in terms of mean sMAPE")
print(mean(rv$sMAPE))
print("performance in terms of mean MASE")
print(mean(rv$MASE))
print("coverage of prediction intervals -- should be close to 95%")
print(mean(rv$InCI)/H)

saveRDS(rv, "tsdl.monthly.homo.rds")

# hetero error
print("monthly tsdl hetero error.....")
s = system.time({
  rv=Rlgt::blgt.multi.forecast(train.data[w.series], future.data[w.series], n.samples=1e4, m = 12, homoscedastic = F)
})
print(s)                         # overall timing info
s[[3]] / length(w.series) # per series time
print("performance in terms of mean sMAPE")
print(mean(rv$sMAPE))
print("performance in terms of mean MASE")
print(mean(rv$MASE))
print("coverage of prediction intervals -- should be close to 95%")
print(mean(rv$InCI)/H)

saveRDS(rv, "tsdl.monthly.hetero.rds")


############################################################################
# quarterly series
H <- 8
train.data = list()
future.data = list()
for (i in 1:length(quarterly.data))
{
  n <- length(quarterly.data[[i]])
  train.data[[i]] = as.numeric(quarterly.data[[i]][1:(n-H)])
  future.data[[i]] = as.numeric(quarterly.data[[i]][(n-H+1):n])
}
w.series = 1:length(quarterly.data)

# homo error
print("quarterly tsdl homo error.....")
s = system.time({
  rv=Rlgt::blgt.multi.forecast(train.data[w.series], future.data[w.series], n.samples=1e4, m = 4, homoscedastic = T)
})
print(s)                         # overall timing info
s[[3]] / length(w.series) # per series time
print("performance in terms of mean sMAPE")
print(mean(rv$sMAPE))
print("performance in terms of mean MASE")
print(mean(rv$MASE))
print("coverage of prediction intervals -- should be close to 95%")
print(mean(rv$InCI)/H)

saveRDS(rv, "tsdl.quarterly.homo.rds")

# hetero error
print("quarterly tsdl hetero error.....")
s = system.time({
  rv=Rlgt::blgt.multi.forecast(train.data[w.series], future.data[w.series], n.samples=1e4, m = 4, homoscedastic = F)
})
print(s)                         # overall timing info
s[[3]] / length(w.series) # per series time
print("performance in terms of mean sMAPE")
print(mean(rv$sMAPE))
print("performance in terms of mean MASE")
print(mean(rv$MASE))
print("coverage of prediction intervals -- should be close to 95%")
print(mean(rv$InCI)/H)

saveRDS(rv, "tsdl.quarterly.hetero.rds")