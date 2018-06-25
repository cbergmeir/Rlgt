library("Mcomp")



sMAPE <- function(forecast, actual ){
  value <- 0
  for (i in 1:length(forecast)){
    value <- value + 2*abs(forecast[i]-actual[i])/(forecast[i]+actual[i])
  } 
  value <- value/length(forecast)*100
  return (value)
}

MASE <- function(forecast, actual ){
  up <- 0
  down <- 0
  for (i in 1:length(forecast)){
    up <- up + abs(forecast[i]-actual[i])
  } 
  up <- up/length(forecast)
  
  for (j in (frequency(forecast)+1):length(forecast)){
    down <- down + abs(actual[j]-actual[j-frequency(forecast)])
  } 
  down <- down/(length(forecast)-frequency(forecast))
  
  return (up/down)
}





#### Testing for sdampen


forval <- list()
for (iter in 646:2829){
  testval[[iter-645]] <-  M3Forecast[["DAMPEN"]][iter]
}

# Data
M3.data <- c(subset(M3,"quarterly"), subset(M3,"monthly"))
error <- 0

testval <- list()
for (iter in 1:length(M3.data)){
  testval[[iter]] <- M3.data[[iter]]$xx
}




## sMape
agg_sMAPE <- 0
corr <- 0
for (k in 1:length(testval)){
  agg_sMAPE <- agg_sMAPE + tryCatch(sMAPE (forval[[k]], testval[[k]]), 
                                    error=function(e) {return (sMAPE (forval[[k]]$mean, testval[[k]])) })
  #agg_sMAPE <- agg_sMAPE + sMAPE (forval[[k]], testval[[k]])
}

agg_sMAPE <- agg_sMAPE/length(testval)


## MASE
agg_MASE <- 0

for (k in 1:length(testval)){
  agg_MASE <- agg_MASE + tryCatch(MASE (forval[[k]], testval[[k]]), 
                                  error=function(e) {return (MASE (forval[[k]]$mean, testval[[k]])) })
}

agg_MASE <- agg_MASE/length(testval)
