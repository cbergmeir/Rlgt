#' @title Rlgt forecast
#' @description  produce forecasts from an \code{\link{rlgtfit}} object
#' @param object rlgtfit object
#' @param xreg input regression matrix
#' @param h forecasting horizon (the default is 10 for annual and 2*periods otherwise)
#' @param level confidence levels for prediction intervals a.k.a. coverage percentiles. Musat be between 0 and 100.
#' @param NUM_OF_TRIALS number of simulations to run. Suggested range is between (1000,5000), but it needs 
#' to be higher for good coverage for very high levels, e.g. 99.8. 
#' @param ... currently not used
#' @return returns a forecast object compatible with the forecast package in R
#' @S3method forecast rlgtfit
#' @method forecast rlgtfit
#' @importFrom forecast forecast 
#' @examples 
#' \dontrun{
#' rlgt_model <- rlgt(lynx, 
#'        control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=2000), verbose=TRUE)
#'
#' # print the model details
#' print(rlgt_model)
#' 
#' # Produce Forecasts for the next 10 years
#' forecast_result <- forecast(rlgt_model, h = 10, level=c(80, 95, 98))
#' 
#' plot(forecast_result,main="Forecasting lynx dataset with LGT model")
#' }
#'
#' @export

 
#object=regModel; xreg=regTrain; level=c(80,95); NUM_OF_TRIALS=2000 
#object=rstanmodel; h = length(actuals); level=c(80,95); NUM_OF_TRIALS=2000
#library(sn)
forecast.rlgtfit <- function(object, 
                             xreg=NULL,
                             h=ifelse(frequency(object$x)>1, 
                                      2*frequency(object$x), 10),
                             level=c(80,95),
                             NUM_OF_TRIALS=2000, ...) {
  
	if (!is.null(xreg) && is.null(dim(xreg))) { # convert non-matrix to matrix
		 xreg <- as.matrix(xreg)
	}												 
  # check if you have non-null and names matched xreg...
  has.regression <- !is.null(xreg)
  # this is different from object$use.regression.  
  # use.regression makes final decision whether to use use.regression
  use.regression <- FALSE
  if (has.regression) {
    # regression components detected
    if (object$use.regression) {
      use.regression <- TRUE
      # override horizon input h here
      h <- nrow(xreg)
    } else {
      message("Current model do not support regression. Regression variables will be ignored.")  #actually, all models currently support regression, but maybe leave for future
    }
  } else {
    if (object$use.regression) {
      stop("Model expects a regression component.")
    }
  } 
  
  if (any(level>100) || any(level<0)) {
    message("Warning: levels must be between 0 and 100. Assuming defaults.")
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
  
  seasonality <- object$seasonality
	seasonality_int=as.integer(seasonality)
	fractSeasonality=seasonality-seasonality_int
	
  seasonality2 <- object$seasonality2
	seasonality2_int=as.integer(seasonality2)
	fractSeasonality2=seasonality2-seasonality2_int
	
	#useGeneralizedSeasonality <- object$useGeneralizedSeasonality  this actually not needed, we check if is.null(powSeasonS)
	MAX_VAL=object$control$MAX_VAL
	MIN_VAL=object$control$MIN_VAL
  
  out <- list(model=object,x=object$x)
  #' @importFrom stats tsp
  
  # start.f is the next(first) forecast period
	tspx <- tsp(out$x)
	if (inherits(out$x,'msts') && !is.null(tspx)) {
			start.f <- tspx[2] + 1/frequency(out$x)  
	} else if (inherits(out$x,'ts') && !is.null(tspx)) {
		#' @importFrom stats tsp
		start.f <- tspx[2] + 1/frequency(out$x)
	} else {#this can be seasonal, but numeric, with seasonality(ies) specified in control
		start.f <- length(out$x)+1
	}
	
  
  # extracting all of the params
  nu <- object$params[["nu"]]
  lastB <- object$params[["lastB"]]
  lastSmoothedInnovSize <- object$params[["lastSmoothedInnovSize"]]
  powx <- object$params[["powx"]]
  
  if(use.regression){
    regCoef <- object$params[["regCoef"]]
		regOffset <- object$params[["regOffset"]]
    if(ncol(xreg) != ncol(regCoef)){
      stop("Error: Number of regression coefficients does not match matrix supplied.")
    }
  }
  
  #these initializations are important, do not remove. (at least some of them :-)
  nuS=Inf; bSmS=0; bS=0; locTrendFractS=0; #t=1; irun=1
  
  if (seasonality>1) {
    s  <- object$params[["s"]]
    sS <- rep(1,seasonality_int+h+1)
  }
  if (seasonality2>1) {
		if (seasonality>seasonality2) {
			stop("seasonality has to be smaller than seasonality2")
			return (NULL)  #for god measure
 		}
    s2  <- object$params[["s2"]]
    sS2 <-rep(1,seasonality2_int+h+1)
		
		if (object$levelMethodId==1) {
			recentVals=rep(0,seasonality_int+h)
			recentVals[1:seasonality_int]=object$x[(length(object$x)-seasonality_int+1):length(object$x)]	
		} else if (object$levelMethodId==2) {
			recentVals=rep(0,seasonality2_int+h)
			recentVals[1:seasonality2_int]=object$x[(length(object$x)-seasonality2_int+1):length(object$x)]	
		}
  }
  
  # Initialise a matrix which contains the last level value
  yf <- matrix(0,nrow=NUM_OF_TRIALS, ncol=h)
  
  # For each forecasting trial
  for (irun in 1:NUM_OF_TRIALS) {
    # Obtain the relevant parameters & use bootstrap sampling
    indx <- sample(nrow(object$params[["l"]]),1)
    prevLevel <- object$params[["lastLevel"]][indx]
    levSmS <- object$params[["levSm"]][indx]
    if (use.regression) {
      regCoefS <- regCoef[indx,]
			regOffsetS <- regOffset[indx]
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
		offsetsigmaS <- object$params[["offsetSigma"]][indx]
    if (!is.null(lastSmoothedInnovSize)) {
      innovSize <- lastSmoothedInnovSize[indx] # we are not updating innovSize in simulation 
      #innovSmS <- object$params[["innovSm"]][indx]
    } else if (!is.null(powx)) {
      powxS <- powx[indx]      
    }
    
    #t=1; irun=1
		if (use.regression) {
			r=xreg %*% regCoefS + regOffsetS  #r = (xreg[t,:] * regCoef + regOffset)*USE_REGRESSION;	
		} else {
			r=rep(0,h)
		}
    if (seasonality > 1) { #seasonal
      powSeasonS <- object$params[["powSeason"]][indx]
			if (fractSeasonality>0) {
				sS[1:(seasonality_int+1)]=s[indx,(ncol(s)-seasonality_int):ncol(s)]
			} else {
				sS[1:seasonality_int]=s[indx,(ncol(s)-seasonality_int):(ncol(s)-1)]  #last element is empty	
			}
      if (seasonality2>1) {
				if (fractSeasonality2>0) {
					sS2[1:(seasonality2_int+1)]=s2[indx,(ncol(s2)-seasonality2_int):ncol(s2)]
				} else {
					sS2[1:seasonality2_int]=s2[indx,(ncol(s2)-seasonality2_int+1):ncol(s2)]	
				}
				powSeasonS2 <- object$params[["powSeason2"]][indx]
				if(use.regression && object$levelMethodId>0){
					regFittedS<-object$params[["r"]][indx]
					if (object$levelMethodId==1) {
						recentVals[1:seasonality_int]=recentVals[1:seasonality_int]-regFittedS[(length(regFittedS)-seasonality_int+1):length(regFittedS)]
					} else {
						recentVals[1:seasonality2_int]=recentVals[1:seasonality2_int]-regFittedS[(length(regFittedS)-seasonality2_int+1):length(regFittedS)]
					}
				}
      }
			
      for (t in 1:h) {
				if (seasonality2>1) {
					if (is.null(powSeasonS)) {
						seasonA <- sS[t]*sS2[t]
					} else {
						seasonA <- sS[t]*abs(prevLevel)^powSeasonS + sS2[t]*abs(prevLevel)^powSeasonS2
					}
				} else if (seasonality > 1) {
					if (is.null(powSeasonS)) {
						seasonA <- sS[t]
					} else {
						seasonA <- sS[t]*abs(prevLevel)^powSeasonS 
					}	
				}

        ## From eq.6.a
        if (is.null(powSeasonS)) {
          expVal <-(prevLevel + coefTrendS*abs(prevLevel)^powTrendS)* seasonA + r[t];	
        } else {
          expVal <- prevLevel + coefTrendS*abs(prevLevel)^powTrendS + seasonA + r[t];
        }
        
        if (!is.null(powx)) {
          omega <- sigmaS*(abs(expVal))^powxS+offsetsigmaS
        } else if (!is.null(lastSmoothedInnovSize)) {
          omega <- sigmaS*innovSize + offsetsigmaS
        } else {
          omega <- sigmaS
        }
        
        # Generate the t-dist error
        error <- rst(n=1, xi=0 ,omega=omega, alpha=0, nu=nuS)
        
        # Fill in the matrix of predicted value
        yf[irun,t] <- min(MAX_VAL,max(MIN_VAL,expVal+error))
        
        # find the currLevel
        if (seasonality2>1) {
					if (is.null(powSeasonS)) {
						if (object$levelMethodId==0) {
							newLevelP=(yf[irun,t]-r[t])/seasonA
						} else if (object$levelMethodId==1) {
							recentVals[seasonality_int+t]=yf[irun,t]-r[t]
							newLevelP=mean(recentVals[(t+1):(t+seasonality_int)])/sS2[t]
						} else if (object$levelMethodId==2) {
							recentVals[seasonality2_int+t]=yf[irun,t]-r[t]
							newLevelP=mean(recentVals[(t+1):(t+seasonality2_int)])
						}	
					} else { #generalized seasonality
						if (object$levelMethodId==0) {
							newLevelP=yf[irun,t]-r[t]-seasonA
						} else if (object$levelMethodId==1) {
							recentVals[seasonality_int+t]=yf[irun,t]-r[t]
							newLevelP=mean(recentVals[(t+1):(t+seasonality_int)]) - sS2[t]*abs(prevLevel)^powSeasonS2
						} else if (object$levelMethodId==2) {
							recentVals[seasonality2_int+t]=yf[irun,t]-r[t]
							newLevelP=mean(recentVals[(t+1):(t+seasonality2_int)])
						}							
					}
					currLevel=max(MIN_VAL,levSmS*newLevelP + (1-levSmS)*prevLevel) ;
				} else if (is.null(powSeasonS)){
					currLevel <- max(MIN_VAL,levSmS*(yf[irun,t]-r[t])/seasonA + (1-levSmS)*prevLevel) ;
				} else {
					currLevel <- max(MIN_VAL,levSmS*(yf[irun,t]-r[t]-seasonA) + (1-levSmS)*prevLevel) ;
				}
    
        if (currLevel>MIN_VAL) {
          prevLevel <- currLevel
        } 
				if (fractSeasonality>0) {
					sS[t+seasonality_int+1] <- sS[t]
					sS[t+seasonality_int]=fractSeasonality*sS[t+seasonality_int]+(1-fractSeasonality)*sS[t]
				} else {
					sS[t+seasonality_int] <- sS[t];	
				}
        
        if (seasonality2>1) {
					if (fractSeasonality2>0) {
						sS2[t+seasonality2_int+1] <- sS2[t]
						sS2[t+seasonality2_int]=fractSeasonality2*sS2[t+seasonality2_int]+(1-fractSeasonality2)*sS2[t]
					} else {
						sS2[t+seasonality2_int] <- sS2[t];	
					}
        }
      }	#through horizons
      # yf[irun,]
    } else { #nonseasonal
      for (t in 1:h) {
        expVal <- prevLevel + coefTrendS*(abs(prevLevel)) ^ powTrendS + locTrendFractS * bS + r[t]
        if (!is.null(powx)) {
					omega <- sigmaS*(abs(expVal))^powxS+offsetsigmaS
        } else if (!is.null(lastSmoothedInnovSize)) {
          omega <- sigmaS * innovSize + offsetsigmaS
        } else {
          omega <- sigmaS
        }
        error <- rst(n=1, xi=0, omega=omega, alpha=0, nu=nuS)
        
        yf[irun,t] <- min(MAX_VAL,max(MIN_VAL, expVal + error))

        ## update level equation
        currLevel <- max(MIN_VAL, levSmS * (yf[irun,t]-r[t]) + (1-levSmS)*prevLevel) ;  #l[t] = levSm*(y[t]-r) + (1-levSm)*l[t-1]; 
				
        ## update trend equations
        if (currLevel>MIN_VAL){
          bS <- bSmS * (currLevel - prevLevel) + (1 - bSmS) * locTrendFractS * bS 
          prevLevel <- currLevel
        } 
        
      } # through horizons
    } # end of seasonal and non-seasonal split
    # extract the with regression components
  } #through trials (simulations)
  
  out$yf <- yf

	#' @importFrom stats quantile
	avgYfs <- apply(yf,2,quantile,probs=quantiles)
	#out$mean <- apply(yf, 2, mean)
	out$mean <- avgYfs[indexOfMedian,]  # Median is safer. We want to be compatible with Forecast package, but there Point Forecast==mean
	if (indexOfMedian > 1) {
	  out$lower <- t(avgYfs[1:(indexOfMedian-1),])
		out$upper <- t(avgYfs[(indexOfMedian+1):(indexOfMedian+length(upperPercentiles)),])
	}
	
	if (inherits(out$x,'msts') && !is.null(tspx)) {
		#' @importFrom forecast msts
		out$mean <- msts(out$mean, seasonal.periods=attributes(out$x)$msts, ts.frequency=tspx[3], start=start.f)  #it is better to get it from the input series and not rely on seasonality(2)  
		if (indexOfMedian > 1) {
			out$lower <- msts(out$lower, seasonal.periods=attributes(out$x)$msts, ts.frequency=tspx[3], start=start.f)
			out$upper <- msts(out$upper, seasonal.periods=attributes(out$x)$msts, ts.frequency=tspx[3], start=start.f)
		}	 
	} else if (inherits(out$x,'ts') && !is.null(tspx)) {
		#' @importFrom stats ts
		out$mean <- ts(out$mean, frequency=seasonality, start=start.f) 
		if (indexOfMedian > 1) {
			out$lower <- ts(out$lower, frequency=seasonality, start=start.f)
			out$upper <- ts(out$upper, frequency=seasonality, start=start.f)
		}		
	} 
	
  out$level <- level
	#str(out, max.level=1)
    
  class(out) <- "forecast"
  out
}
