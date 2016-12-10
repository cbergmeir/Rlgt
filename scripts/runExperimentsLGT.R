library("forecast")
library("Mcomp")

library("Rlgt")

#M3.data <- append(subset(M3,"yearly"), append(subset(M3,"quarterly"), subset(M3,"monthly")))
M3.data <- subset(M3,"yearly")
#M3.data <- subset(M3,"monthly")
nseries <- length(M3.data)

#M3.data[[2]]
#M3[["N0001"]]
#str(M3)

mySgeApply <- lapply

#library(parallel)
#mySgeApply <- function(...) { mclapply(..., mc.cores=4)}

set.seed(8)

stanModLGT <- init.lgt()

#nseries <- 5

#curr_series <- 1

#options(error=recover)

#1:nseries
#
forecasts <- mySgeApply(1:nseries, function(curr_series) {
      
      cat(curr_series, "\n")
      
      sizeTestSet <- length(M3.data[[curr_series]]$xx)
      
      data.train <- M3.data[[curr_series]]$x
      
      mod <- list()
      forecasts <- list()
      
#benchmarks
#      mod[["etsD"]] <- ets(data.train, model="AAN", damped=TRUE)
#      forecasts[["etsD"]] <- forecast(mod[["etsD"]], PI=FALSE, h=sizeTestSet)$mean
#      
#      mod[["ets"]] <- ets(data.train, model="AAN")
#      forecasts[["ets"]] <- forecast(mod[["ets"]], PI=FALSE, h=sizeTestSet)$mean
#
#      mod[["etsDcmaes"]] <- ets(data.train, model="AAN", damped=TRUE, solver="malschains_c", 
#          control=malschains.control(ls="cmaes", lsOnly=TRUE))
#      forecasts[["etsDcmaes"]] <- forecast(mod[["etsDcmaes"]], PI=FALSE, h=sizeTestSet)$mean
#      
#      mod[["etsDMalsCh"]] <- ets(data.train, model="AAN", damped=TRUE, solver="malschains_c", 
#          control=malschains.control(ls="cmaes"), maxit=10000)
#      forecasts[["etsDMalsCh"]] <- forecast(mod[["etsDMalsCh"]], PI=FALSE, h=sizeTestSet)$mean

#      set.seed(curr_series)
#      mod[["baggedETS"]] <- baggedETS(data.train)
#      forecasts[["baggedETS"]] <- forecast(mod[["baggedETS"]], h=sizeTestSet)$mean

      set.seed(curr_series)
      mod[["ets"]] <- ets(data.train)
      forecasts[["ets"]] <- forecast(mod[["ets"]], h=sizeTestSet)$mean
      
      #--------------------------------
      #Fit LGT model
      
      mod[["lgt"]] <- fit.lgt(data.train, h = sizeTestSet, stanModel=stanModLGT, ncores=4)
      forecasts[["lgt"]] <- forecast(mod[["lgt"]], h = sizeTestSet)
      
#      
#      set.seed(curr_series)
#      forecasts[["baggedETSold"]] <- forecast:::forecast.baggedETSold(data.train, h=sizeTestSet)$mean
#      
      forecasts
      
    })

print("finished")

#forecast.lgt <- function(x, h) {
#
#  yy=as.numeric(M3[[idr]]$xx)
#  ymax=max(c(y,yy),max(avgYfs[3,]))
#  ymin=min(min(c(y,yy)),min(avgYfs[1,]))
#  plot(c(y,yy), main=series,  col=c(rep(1,n),rep(2,maxPredictionHorizon)), ylab='',ylim=c(ymin,ymax) )
#  for (irun in 1:NUM_OF_TRIALS_TO_SHOW) 
#    lines((n+1):(n+maxPredictionHorizon),yf[irun,], col='gray')
#  lines((n+1):(n+maxPredictionHorizon),avgYfs[1,], col='pink',lwd=2)
#  lines((n+1):(n+maxPredictionHorizon),avgYfs[2,], col='blue',lwd=2)
#  lines((n+1):(n+maxPredictionHorizon),avgYfs[3,], col='pink',lwd=2)
#  lines(c(y,yy), main=series,  col=c(rep(1,n),rep(2,maxPredictionHorizon)),type='b')
#  
#}
#
#
#
#
#
#
#
#
#


#mod[["etsDcmaes"]]
#str(forecasts)
#forecasts

#now run evalExperiments