library("forecast")
library("Mcomp")

library("Rlgt")
library("RlgtLik")

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

#set.seed(8)
set.seed(5)
#stanModLGT <- init.lgt()

#nseries <- 5

curr_series <- 250

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
      mod[["etsAAN"]] <- ets(data.train, model="AAN")
      forecasts[["etsAAN"]] <- forecast(mod[["etsAAN"]], PI=FALSE, h=sizeTestSet)$mean

      mod[["ets"]] <- ets(data.train)
      forecasts[["ets"]] <- forecast(mod[["ets"]], PI=FALSE, h=sizeTestSet)$mean

      mod[["etsLGT"]] <- etsLGT(data.train, bounds="usual")
      forecasts[["etsLGT"]] <- forecast(mod[["etsLGT"]], PI=FALSE, h=sizeTestSet, simulate=TRUE)$mean
      
      mod[["etsLGTcmaes"]] <- etsLGT(data.train, , bounds="usual", solver="malschains_c", 
          control=malschains.control(ls="cmaes", lsOnly=TRUE))
      forecasts[["etsLGTcmaes"]] <- forecast(mod[["etsLGTcmaes"]], PI=FALSE, h=sizeTestSet, simulate=TRUE)$mean
     
      
      
#      mod[["etsDMalsCh"]] <- etsLGT(data.train, model="AAN", damped=TRUE, solver="malschains_c", 
#          control=malschains.control(popsize = 5000, ls="cmaes"), maxit=50000)
#      forecasts[["etsDMalsCh"]] <- forecast(mod[["etsDMalsCh"]], PI=FALSE, h=sizeTestSet, simulate=TRUE)$mean

#      set.seed(curr_series)
#      mod[["baggedETS"]] <- baggedETS(data.train)
#      forecasts[["baggedETS"]] <- forecast(mod[["baggedETS"]], h=sizeTestSet)$mean
#
#      set.seed(curr_series)
#      mod[["ets"]] <- ets(data.train)
#      forecasts[["ets"]] <- forecast(mod[["ets"]], h=sizeTestSet)$mean
#      
      #--------------------------------
      #Fit LGT model
      
#      mod[["lgt"]] <- fit.lgt(data.train, stanModel=stanModLGT, ncores=4)
#      forecasts[["lgt"]] <- forecast(mod[["lgt"]], h = sizeTestSet)
#      
#      
#      set.seed(curr_series)
#      forecasts[["baggedETSold"]] <- forecast:::forecast.baggedETSold(data.train, h=sizeTestSet)$mean
#      
      forecasts
      
    })

print("finished")










data.test <- M3.data[[curr_series]]$xx

#names(mod[["etsDcmaes"]])
#
#mod[["etsDcmaes"]]$par
#
#mod[["lgt"]]
#
#par(mfrow=c(2,1))
#plot(forecast(mod[["etsDcmaes"]], h=6))
#lines(data.test)
#
#plot(forecast(mod[["lgt"]], h=6))
#lines(data.test)
##now run evalExperiments

#names(res.claw)

stanVec <- mod[["lgt"]]$paramMeans

mod[["etsD"]]$par
mod[["lgt"]]$paramMeans

oldPar <- mod[["etsDcmaes"]]$par
oldState <- mod[["etsDcmaes"]]$state

oldParams <- mod[["lgt"]]$params

#obj$state <- oldState

newPar <- c(alpha=stanVec$levSm, beta=stanVec$bSm,
    phi=stanVec$locTrendFract, lambda=stanVec$coefTrend,
    rho=stanVec$powTrend, l=stanVec$l[1], b=stanVec$b[1])

newState <- t(rbind(stanVec$l, stanVec$b))

#TODO: why are vectors l and b one value shorter?? this could be a problem
#currently, I just duplicate the first value
newState <- rbind(newState[1,], newState)
colnames(newState) <- c("l", "b")
newState <- ts(newState)

tspx <- tsp(mod[["etsDcmaes"]]$state)
#tspx[1] <- tspx[1]+1 
tsp(newState) <- tspx



mod[["etsDcmaes"]]$state <- newState
mod[["etsDcmaes"]]$initstate <- mod[["etsDcmaes"]]$state[1,]
mod[["etsDcmaes"]]$par <- newPar

mod[["lgt_orig"]] <- mod[["lgt"]] 

mod[["lgt"]]$params <- mod[["lgt"]]$paramMeans
mod[["lgt"]]$params[["l"]] <- t(mod[["lgt"]]$params[["l"]])
mod[["lgt"]]$params[["b"]] <- t(mod[["lgt"]]$params[["b"]])

#pdf(file="/home/bergmeir/20161217_series250_forecasts_LGT.pdf")
par(mfrow=c(2,2))
plot(forecast(mod[["etsDcmaes"]], simulate=TRUE, PI=FALSE, h=sizeTestSet), main="paramMeans, my forecast func")
#plot(forecast(mod[["etsDMalsCh"]], simulate=TRUE, PI=FALSE, h=sizeTestSet), main="1", ylim=c(3000,6000))
lines(data.test)

plot(forecast(mod[["etsD"]], simulate=TRUE, PI=FALSE, h=sizeTestSet), main="my implementation of LGT", ylim=c(3000,6000))
lines(data.test)

plot(forecast(mod[["lgt"]], h=sizeTestSet), main="paramMeans used for forecasting")
lines(data.test)

plot(forecast(mod[["lgt_orig"]], h=sizeTestSet), main="Original LGT")
lines(data.test)
#dev.off()

mod[["etsDMalsCh"]]$par
mod[["lgt_orig"]]$paramMeans

forecast(mod[["etsDcmaes"]], PI=FALSE, h=sizeTestSet, simulate=TRUE, h=6)

forecast(mod[["etsDMalsCh"]], PI=FALSE, h=sizeTestSet, simulate=TRUE)


names(mod[["etsDcmaes"]])



#TODO: forecasts should be more or less the same, but they are not

#obj <- mod[["etsDcmaes"]]
#obj$state[length(obj$x)+1,]
#obj$state
#
#
#obj$x
#
#
