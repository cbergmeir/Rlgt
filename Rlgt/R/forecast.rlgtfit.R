
#' This function produces forecasts from a model
#' 
#' @title produce forecasts
#' @param object lgt object
#' @param h Forecasting horizon (10 for annual and 2*periods otherwise)
#' @param level Confidence levels for prediction intervals a.k.a. coverage percentiles. Beween 0 and 100.
#' @param NUM_OF_TRIALS Number of simulations to run. Suggested rannge (1000,5000), but it may have to be higher for good coverage of very high levels, e.g. 99.8. 
#' @param MIN_VAL Minimum value the forecast can take. Must be positive.
#' @param MAX_VAL Maximum value the forecast can take.
#' @param ... currently not being used
#' @return returns a forecast object compatible with the forecast package
#' @S3method forecast lgt
#' @method forecast lgt
#' @importFrom forecast forecast 
#' @examples 
#' \dontrun{
#' lgt_model <- fit.lgt(lynx, model="LGT", nCores=4, nChains=4,
#' control=lgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=2000), 
#' verbose=TRUE)

#' # print the model details
#' print(lgt_model)
#' 
#' # Produce Forecasts for the next 10 years
#' forecast_result <- forecast(lgt_model, h = 10, level=c(80, 95, 98))
#' 
#' plot(forecast_result,main="Forecasting lynx dataset with LGT model")
#' }
#'
#' \dontrun{demo(exampleScript)}
#' @export

# object=mod[["lgte"]]; level=c(80,95, 98); NUM_OF_TRIALS=2000; MIN_VAL=0.001; MAX_VAL=1e38; h=8
# library(sn)
forecast.lgt <- function(object, 
                         xreg=NULL,
                         h=ifelse(frequency(object$x)>1, 
                                  2*frequency(object$x), 10),
                         level=c(80,95),
                         NUM_OF_TRIALS=2000, 
                         MIN_VAL=0.001, MAX_VAL=1e38, ...) {
  
  # check if you have non-null and names matched xreg...
  # WIP
  if (object$has.regression && is.null(xreg)){
    stop("")
  }
  
  if (any(level>100) || any(level<0)) {
    message("Warning: levels mus be between 0 and 100. Assuming defaults.")
    level=c(80,95)
  }
  
  if (length(level)==1 && level==50) {
    percentiles=level
    indexOfMedian=1
  } else {
    if (50 %in% level) {
      level=level[level!=50]
    } 
    lowerPercentiles=(100-level)/2
    indexOfMedian=length(lowerPercentiles)+1
    upperPercentiles=100-lowerPercentiles
    percentiles=c(lowerPercentiles,50,upperPercentiles) #this follows convention of forecast package where forecast$lower are from higher to lower
  }
  
  quantiles=percentiles/100.
  
  SEASONALITY <- object$control$SEASONALITY
  SEASONALITY2 <- object$control$SEASONALITY2
  
  out <- list(model=object,x=object$x)
  #' @importFrom stats tsp
  #' @importFrom stats tsp<-
  tspx <- tsp(out$x)
  
  # start.f is the next(first) forecast period
  if (!is.null(tspx)) {
    start.f <- tspx[2] + 1/frequency(out$x)
  } else {
    start.f <- length(out$x)+1
  }
  
  # extracting all of the params
  nu <- object$params[["nu"]]
  lastB <- object$params[["lastB"]]
  lastSmoothedInnovSize <- object$params[["lastSmoothedInnovSize"]]
  powx <- object$params[["powx"]]
  
  #these initializations are important, do not remove. 
  nuS=Inf; bSmS=0; bS=0; locTrendFractS=0; #t=1; irun=1
  
  if (SEASONALITY>1) {
    s  <- object$params[["s"]]
    sS <- rep(1,SEASONALITY+h)
  }
  if (SEASONALITY2>1) {
    s2  <- object$params[["s2"]]
    sS2 <-rep(1,SEASONALITY2+h)
  }
  
  # Initialise a matrix which contains the last level value
  yf <- matrix(0,nrow=NUM_OF_TRIALS, ncol=h)
  
  # For each forecasting trial
  for (irun in 1:NUM_OF_TRIALS) {
    # Obtain the relevant parameters & use bootstrap sampling
    indx <- sample(nrow(object$params[["l"]]),1)
    prevLevel <- object$params[["lastLevel"]][indx]
    levSmS <- object$params[["levSm"]][indx]
    if (object$has.regression) {
      regCoef <- object$params[["regCoef"]][indx]
    }
    
    if (!is.null(nu)) {
      nuS <- nu[indx]
    } 
    
    powTrendS  <- object$params[["powTrend"]][indx]
    coefTrendS <- object$params[["coefTrend"]][indx]
    
    if (!is.null(lastB)) {
      bS <- lastB[indx]
      bSmS <- object$params[["bSm"]][indx]
      locTrendFractS <- object$params[["locTrendFract"]][indx]
    }
    
    sigmaS <- object$params[["sigma"]][indx]
    if (!is.null(lastSmoothedInnovSize)) {
      innovSize <- lastSmoothedInnovSize[indx]
      innovSmS <- object$params[["innovSm"]][indx]
      offsetsigmaS <- object$params[["offsetSigma"]][indx]
    } else if (!is.null(powx)) {
      powxS <- powx[indx]
      offsetsigmaS <- object$params[["offsetSigma"]][indx]
    }
    
    #t=1; irun=1
    if (SEASONALITY > 1) { #seasonal
      powSeasonS <- object$params[["powSeason"]][indx]
      
      sS[1:SEASONALITY]=s[indx,(ncol(s)-SEASONALITY+1):ncol(s)]
      if (SEASONALITY2>1) {
        sS2[1:SEASONALITY2]=s2[indx,(ncol(s2)-SEASONALITY2+1):ncol(s2)]
      }
      for (t in 1:h) {
        seasonA=sS[t]
        if (SEASONALITY2>1) {
          seasonA <- seasonA * sS2[t]
        }
        ## From eq.6.a
        if (is.null(powSeasonS)) {
          expVal <- (prevLevel+ coefTrendS * abs(prevLevel) ^ powTrendS) * seasonA;	
        } else {
          expVal <- prevLevel+ coefTrendS * abs(prevLevel) ^ powTrendS + seasonA *
            abs(prevLevel) ^ powSeasonS;
        }
        
        if (!is.null(powx)) {
          omega <- sigmaS*(abs(prevLevel))^powxS+offsetsigmaS
        } else if (!is.null(lastSmoothedInnovSize)) {
          omega <- sigmaS*innovSize + offsetsigmaS
        } else {
          omega <- sigmaS
        }
        
        # Generate the t-dist error
        error <- rst(n=1, xi=0 ,omega=omega, alpha=0, nu=nuS)
        
        # Fill in the matrix of predicted value
        yf[irun,t]=min(MAX_VAL,max(MIN_VAL,expVal+error))
        
        # find the currLevel
        if (inherits(object$model, "RlgtStanModelSGT2")) {
          currLevel=max(MIN_VAL,levSmS*yf[irun,t]/seasonA + (1-levSmS)*expVal/seasonA) ;
        }
        else if (is.null(powSeasonS)){
          currLevel=max(MIN_VAL,levSmS*yf[irun,t]/seasonA + (1-levSmS)*prevLevel) ;
        } else {#gSTG model
          currLevel=max(MIN_VAL,levSmS*(yf[irun,t]-seasonA*abs(prevLevel)^powSeasonS) + (1-levSmS)*prevLevel) ;
        }     
        if (currLevel>MIN_VAL) {
          prevLevel=currLevel
        } 
        sS[t+SEASONALITY] <- sS[t];
        if (SEASONALITY2>1) {
          sS2[t+SEASONALITY2] <- sS2[t];
        }
      }	#through horizons
      # yf[irun,]
    } else { #nonseasonal
      for (t in 1:h) {
        expVal <- prevLevel + coefTrendS*(abs(prevLevel)) ^ powTrendS + locTrendFractS * bS
        if (!is.null(powx)) {
          omega <- sigmaS * (abs(prevLevel)) ^ powxS + offsetsigmaS
        } else if (!is.null(lastSmoothedInnovSize)) {
          omega <- sigmaS * innovSize + offsetsigmaS
        } else {
          omega <- sigmaS
        }
        error <- rst(n=1, xi=0, omega=omega, alpha=0, nu=nuS)
        
        yf[irun,t] <- min(MAX_VAL,max(MIN_VAL, expVal + error))
        
        ## update level equation
        if (inherits(object$model, "RlgtStanModelLGT2")) {
          currLevel <- max(MIN_VAL, levSmS * yf[irun,t] + (1-levSmS) * expVal) ;
        }
        else {
          currLevel <- max(MIN_VAL, levSmS * yf[irun,t] + (1-levSmS)*prevLevel) ;
        }
        
        ## update trend equations
        if (currLevel>MIN_VAL & !inherits(object$model, "RlgtStanModelLGT2")) {
          # but bSmS and bS may be==0 so then noop
          prevLevel <- currLevel
          bS <- bSmS * (currLevel - prevLevel) + (1 - bSmS) * bS 
        } 
        
        else if (currLevel>MIN_VAL){
          # but bSmS and bS may be==0 so then noop
          bS <- bSmS * (currLevel - prevLevel) + (1 - bSmS) * locTrendFractS * bS 
          prevLevel <- currLevel
        } 
        
      } #through horizons
      # yf[irun,]
    }
  } #through trials (simulations)
  
  out$yf <- yf
  
  #' @importFrom stats ts
  out$mean <- ts(apply(yf, 2, mean), frequency=max(frequency(out$x),SEASONALITY), 
                 start=start.f ) #still not right if there are two seasonalities
  
  #' @importFrom stats quantile
  avgYfs=apply(yf,2,quantile,probs=quantiles)
  
  out$median <- ts(avgYfs[indexOfMedian,])
  if (indexOfMedian > 1) {
    out$lower <- ts(t(avgYfs[1:(indexOfMedian-1),]))
    out$upper <- ts(t(avgYfs[(indexOfMedian+1):(indexOfMedian+length(upperPercentiles)),]))	
  }
  
  out$level <- level
  
  # Assign tsp attributes for consistency
  if (SEASONALITY2<=1) {
    tsp(out$median) <- tsp (out$lower) <- tsp(out$upper) <- tsp(out$mean)
  }
  # for dual seasonality this fails with Error in `tsp<-`(`*tmp*`, value = c(653, 700, 1)) : 
  #invalid time series parameters specified
  
  class(out) <- "forecast"
  out
}
