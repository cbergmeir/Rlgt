
#' Initialize a stan model that uses the (non-seasonal) LGT
#' 
#' @title Initialize a non-seasonal LGT stan model
#' @param modelType type of the forecasting model selected
#' @returnType RlgtStanModelSGT
#' @return SkeletonModel
#' 
#' @importFrom rstan stan_model
#' @export
initModel <- function(modelType = NULL){
  
  if(is.null(modelType)) {
    print("No model type was provided, generating an LGT model.")
    modelType <- "LGT"
  }
  
  model <- list()
	
  if(modelType=="LGT") {
    #Non-Seasonal Local Global Trend model
		model[["parameters"]] <- c("l", "b", "nu", "sigma", "levSm",  "bSm", 
				"powx", "coefTrend",  "powTrend", "offsetSigma", "locTrendFract")
    
    model[["model"]] <- stanmodels$lgt
    class(model) <- c("RlgtStanModelLGT")
    
  } else if (modelType=="LGTe") {
    #Non-Seasonal Local Global Trend model with smoothed error size
    model[["parameters"]] <- c("l", "b", "smoothedInnovSize", 
				"coefTrend",  "powTrend", "sigma", "offsetSigma",
				"bSm", "locTrendFract", "levSm", "innovSm", "nu")
    
    model[["model"]] <- stanmodels$LGTe
    class(model) <- c("RlgtStanModelLGTe")
  } else if(modelType=="SGT") {
		#Seasonal Global Trend model
		model[["parameters"]] <- c("l", "s", "sSm","nu", "sigma", "levSm", 
				"powx", "coefTrend", "powTrend", "offsetSigma")
		
		model[["model"]] <- stanmodels$SGT
		class(model) <- c("RlgtStanModelSGT")
		
	} else if (modelType=="SGTe") {
		#Seasonal Global Trend model with smoothed error size 
		model[["parameters"]] <- c("l", "s", "smoothedInnovSize", 
				"coefTrend", "powTrend", "sigma", "offsetSigma",
				"levSm", "sSm", "innovSm", "nu")
		
		model[["model"]] <- stanmodels$SGTe
		class(model) <- c("RlgtStanModelSGTe")
	}	else if (modelType=="Trend") {
		#trend-only, Gaussian error, homoscedastic version of LGT
		model[["parameters"]] <- c("l", "b", "sigma", "levSm",  "bSm", "coefTrend", 
				"powTrend", "locTrendFract")
		
		model[["model"]] <- stanmodels$trend
		class(model) <- c("RlgtStanModelTrend")
	} 
  
  class(model) <- c("RlgtStanModel", class(model))
  model  
  
} 


#' @title lgt class
#' @description a constructor function for the "lgt" class
#' @param y the time series data
#' @param lgtmodel type of lgtmodel selected
#' @param params list of parameters
#' @param paramMean mean of each parameter
#' @param seasonality number of seasons, 1 for annual
#' @param samples stanfit object representing the MCMC samples
#' @return lgt instance

lgt <- function(y,lgtmodel,params, paramMean, seasonality, samples) {
  # we can add our own integrity checks

  value <- list(x = y, model = lgtmodel, params = params, paramMeans=paramMean, SEASONALITY=seasonality, samples=samples)
  
  # class can be set using class() or attr() function
  attr(value, "class") <- "lgt"
  value
}

#### This function is created mainly to provide a default value for each parameter

#'
#' @title Sets and initializes the main parameters of the algorithm
#' @description This is a function that initializes and sets the parameters of the algorithm. 
#' It generates a list of parameters, to be used with the \code{\link{fit.lgt}} function. 
#' @param MAX_RHAT_ALLOWED Maximum average Rhat that suggests a good fit, see Stan's manual. Suggested range(1.005,1.02), see also MAX_NUM_OF_REPEATS description below.
#' @param NUM_OF_ITER Number of iterations for each chain. Suggested range(1000,5000). Generally, the longer the series, the smaller the vallue will do. 
#' See also MAX_NUM_OF_REPEATS description below.
#' @param MAX_NUM_OF_REPEATS Maximum number of the sampling procedure repeats if the fit is unsatisfactorily (avgRHat>MAX_RHAT_ALLOWED).
#' Each round doubles the number of iterations. Suggested range(2,4)
#' @param CAUCHY_SD_DIV For parameters with non-obvious range Cauchy distribution is used. The error size of this distribution 
#' is calculated by dividing max value of the time series by this constant. Suggested range(100,300)
#' @param MIN_SIGMA Minimum size of the fitted sigma, applied for numerical stability. Must bve positive. 
#' @param MIN_NU Minimum degrees of freedom of the Student's distribution, that is used in most models. Suggested range(1.2, 5)
#' @param MAX_NU Maximum degrees of freedom of the Student's distribution. Suggested range(15,30) 
#' @param MIN_POW_TREND Minimum value of power of trend coefficient. Suggested range(-1,0) 
#' @param MAX_POW_TREND Maximum value of power of trend coefficient. It should stay 1 to allow the model to approach exponential growth when needed.
#' @param POW_TREND_ALPHA Alpha parameter of Beta distribution that is the prior of the power coefficient in the formula of trend parameter.
#' To make the forecast more curved, make it larger. Suggested range(1,6)
#' @param	POW_TREND_BETA Beta parameter of Beta distribution that is the prior of the power of trend parameter. 1 by default, see also above.
#' @param	POW_SIGMA_ALPHA Alpha parameter of Beta distribution that is the prior of the power coefficient in the formula of the error size. 1 by default, see also below.
#' @param	POW_SIGMA_BETA Beta parameter of Beta distribution that is the prior of the power coefficient in the formula of the error size.
#' If the powSigma fitted is considered too often too high (i.e.> 0.7) you can attempt to tame it down by increasing POW_SIGMA_BETA.  Suggested range(1,4).
#' ADAPT_DELTA Target Metropolis acceptance rate. See Stan manual. Suggested range (0.8-0.97). 
#' MAX_TREE_DEPTH NUTS maximum tree depth. See Stan manual. Suggested range (10-12).
#' @param	SEASONALITY E.g. 12 for monthly seasonality. 1 for non-seasonal models
#' @param	SKEW Skew of error distribution used by manually-skewed models. 0 be default. 
#' @param MAX_TREE_DEPTH Description
#' @param ADAPT_DELTA Description
#' Setting it negative makes negative innovations having smaller impact on the fitting than the positive ones,
#' which would have the effect of making a model "more optimistic". Suggested range (-0.5, 0.5).
#' @returnType 
#' @return 
#' @export
lgt.control <- function(
    MAX_RHAT_ALLOWED=1.005, 
    NUM_OF_ITER=2500,
    MAX_NUM_OF_REPEATS=3,
    CAUCHY_SD_DIV=200,
    MIN_SIGMA=0.001, 
    MIN_NU=2,
    MAX_NU=20,
    MIN_POW_TREND=-0.5,
    MAX_POW_TREND=1,
		POW_TREND_ALPHA=1,
		POW_TREND_BETA=1,
		POW_SIGMA_ALPHA=1,  
		POW_SIGMA_BETA=1, 
		ADAPT_DELTA=0.9, 
		MAX_TREE_DEPTH=11,
		SEASONALITY=1,
		SKEW=0 
) {
  
  list(CAUCHY_SD_DIV=CAUCHY_SD_DIV,
      MAX_RHAT_ALLOWED=MAX_RHAT_ALLOWED,
      NUM_OF_ITER=NUM_OF_ITER,
      MAX_NUM_OF_REPEATS=MAX_NUM_OF_REPEATS,
      MIN_SIGMA=MIN_SIGMA,
      MIN_NU=MIN_NU,
      MAX_NU=MAX_NU,
      MIN_POW_TREND=MIN_POW_TREND,
      MAX_POW_TREND=MAX_POW_TREND,
			POW_TREND_ALPHA=POW_TREND_ALPHA,
			POW_TREND_BETA=POW_TREND_BETA,
			POW_SIGMA_ALPHA=POW_SIGMA_ALPHA,
			POW_SIGMA_BETA=POW_SIGMA_BETA,
			ADAPT_DELTA=ADAPT_DELTA,
			MAX_TREE_DEPTH=MAX_TREE_DEPTH,
			SEASONALITY=SEASONALITY,
			SKEW=SKEW
			)
}




#' Runs the model fitting
#' 
#' @title Runs the model fitting
#' @param y the time series
#' @param model a stan model
#' @param control control arguments list
#' @param nChains number of MCMC chains . Must >=1. Perhaps optimal number is 4.
#' @param nCores number of cores to be used. For performance reasons it should be equal to nChains, 
#' but nChains should be smaller or equal to the number of cores on the computer.  
#' @param addJitter adding a bit of jitter is helping Stan in case of some flat series
#' @param verbose print verbose information yes/no
#' @returnType lgt
#' @return lgtModel
#' 
#' @importFrom rstan rstan_options
#' @importFrom rstan sampling
#' @importFrom rstan get_sampler_params
#' @importFrom rstan extract
#' @importMethodsFrom rstan summary
#' @importFrom sn rst
#' @export
fit.lgt <- function(y, model=c("LGT", "SGT", "LGTe", "SGTe", "Trend"), 
	control=lgt.control(), nChains=2, nCores=2, addJitter=TRUE, verbose=FALSE) {
  
  if(!inherits(model, "RlgtStanModel")) {
    model <- initModel(model)
  }
  
  MAX_RHAT_ALLOWED=control$MAX_RHAT_ALLOWED
  NUM_OF_ITER=control$NUM_OF_ITER
  MAX_NUM_OF_REPEATS=control$MAX_NUM_OF_REPEATS
  
  if(nCores>1) {
    rstan_options(auto_write = TRUE)
    options(mc.cores = nCores)    
  }
  
  y.orig <- y
  n=length(y)
  
  #' @importFrom stats rnorm
  if(addJitter) y <- y + rnorm(n,0,sd=abs(min(y))*0.0001)
  
  CauchySd=max(y)/control$CAUCHY_SD_DIV
	
	SEASONALITY=control$SEASONALITY
	#' @importFrom stats frequency
	if (frequency(y)>1) SEASONALITY=frequency(y) #good idea?
	
	data <- list(CAUCHY_SD=CauchySd, SEASONALITY=SEASONALITY,
			MIN_POW_TREND=control$MIN_POW_TREND, 
			MAX_POW_TREND=control$MAX_POW_TREND, 
			MIN_SIGMA=control$MIN_SIGMA,  
			MIN_NU=control$MIN_NU,  
			MAX_NU=control$MAX_NU, 
			POW_SIGMA_ALPHA=control$POW_SIGMA_ALPHA, 
			POW_SIGMA_BETA=control$POW_SIGMA_BETA,
			POW_TREND_ALPHA=control$POW_TREND_ALPHA, 
			POW_TREND_BETA=control$POW_TREND_BETA,
			y=y, N=n, SKEW=control$SKEW) # to be passed on to Stan
	
  ### Repeat until Rhat is in an acceptable range (i.e. convergence is reached)
  avgRHat=1e200; irep=1
  for (irep in 1:MAX_NUM_OF_REPEATS) {
		initializations <- list();
		for (irr in 1:nChains) {
			initializations[[irr]]=list( 
			  ### Initialise seasonality factors
					initSu=rnorm(SEASONALITY,1,0.1) # for non-seasonal models it is not necessary, but makes code simpler and is not a big overhead
			)
		}

		#Double the number of iterations for the next cycle
		numOfIters=NUM_OF_ITER*2^(irep-1)
    samples1=
      sampling(
				control=list(adapt_delta = control$ADAPT_DELTA, max_treedepth=control$MAX_TREE_DEPTH),
        model$model,   
        data=data, 
        init=initializations,
        pars=model$parameters,
        iter=numOfIters,
        chains=nChains,
        cores=nCores,
        open_progress=F,
        refresh = if(verbose) numOfIters/5 else -1)

    ### Get the Rhat values
    ainfo=summary(samples1)
    RHats=ainfo$summary[,10]
    RHats=as.numeric(RHats[is.finite(RHats)])
    currRHat=mean(RHats, na.rm=T)
    if (currRHat<=MAX_RHAT_ALLOWED) {
      samples=samples1
      avgRHat=currRHat
      print(samples)
      print(paste("avgRHat",avgRHat))
      break
    } else {
      if (currRHat<avgRHat) {#in the worst case this is at least once executed, because avgRHat is initialized high
        samples=samples1
        avgRHat=currRHat
        print(samples)
        print(paste("avgRHat",avgRHat))
      } else {
        print ("worse...")
        print(paste("currRHat",currRHat))
      }
      if (MAX_NUM_OF_REPEATS>irep) print (paste("trying to do better..."))
    }
    #str(samples, max.level =4)
  }#repeat if needed
  
  if(verbose) print(summary(do.call(rbind, args = get_sampler_params(samples1, inc_warmup = F)), digits = 2)) #diagnostics including step sizes and tree depths
  
  params <- list()
  paramMeans <- list()
  
  # Extract all of the parameter means
  for(param in model$parameters) {
    params[[param]] <- extract(samples)[[param]]
    
    # Find the mean, but for ETS components average based on each t
    paramMeans[[param]] <- if(param %in% c("l", "b", "s")) apply(params[[param]],2,mean) else mean(params[[param]])
  }
	
	#special processing for vector params, where only the last value counts
  # This includes level, trend, and innov size
  params[["lastLevel"]] <- params[["l"]][,ncol(params[["l"]])]
	paramMeans[["lastLevel"]] <- mean(params[["lastLevel"]])
	
	if ("b" %in% model$parameters) {
		params[["lastB"]] <- params[["b"]][,ncol(params[["b"]])]
		paramMeans[["lastB"]] <- mean(params[["lastB"]])
	}
	if ("smoothedInnovSize" %in% model$parameters) {
		params[["lastSmoothedInnovSize"]] <- params[["smoothedInnovSize"]][,ncol(params[["smoothedInnovSize"]])]
		paramMeans[["lastSmoothedInnovSize"]] <- mean(params[["lastSmoothedInnovSize"]])
	}

	
  out <- lgt(y.orig, model, params, paramMeans, SEASONALITY, samples)

  out

}




#' This function produces forecasts from a model
#' 
#' @title produce forecasts
#' @param object lgt object
#' @param h Forecasting horizon (10 for annual and 2*periods otherwise)
#' @param level Confidence levels for prediction intervals a.k.a. coverage percentiles. Beween 0 and 100.
#' @param NUM_OF_TRIALS Number of simulations to run. Suggested rannge (1000,5000), but it may have to be higher for good coverage of very high levels, e.g. 99.8. 
#' @param MIN_VAL Minimum value the forecast can take. Must be positive.
#' @param MAX_VAL Maximum value the forecast can take.
#' @param ... description
#' @returnType forecast
#' @return returns a forecast object compatible with the forecast package
#' @S3method forecast lgt
#' @method forecast lgt
#' @importFrom forecast forecast 
#' @author bergmeir
#' @export

# object=mod[["lgte"]]; level=c(80,95, 98); NUM_OF_TRIALS=2000; MIN_VAL=0.001; MAX_VAL=1e38; h=8
# library(sn)
forecast.lgt <- function(object, h=ifelse(frequency(object$x)>1, 2*frequency(object$x), 10),
		level=c(80,95),
		NUM_OF_TRIALS=2000, 
    MIN_VAL=0.001, MAX_VAL=1e38, ...) {

  #object <- mod[["lgt"]]
	
	if (any(level>100) || any(level<0)) {
		print(paste("Warning: levels mus be between 0 and 100. Assuming defaults."))
		level=c(80,95)
	}
	
	if(length(level)==1 && level==50) {
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
	
	SEASONALITY <- object$SEASONALITY 
  
  out <- list(model=object,x=object$x)
  #' @importFrom stats tsp
  #' @importFrom stats tsp<-
  tspx <- tsp(out$x)
  
  # start.f is the next(first) forecast period
  if(!is.null(tspx)){
		start.f <- tspx[2] + 1/frequency(out$x)
	} else {
		start.f <- length(out$x)+1
	}
    
  # extracting all of the params
	nu=object$params[["nu"]]
	lastB <- object$params[["lastB"]]
	lastSmoothedInnovSize <- object$params[["lastSmoothedInnovSize"]]
	powx <- object$params[["powx"]]
	s <- object$params[["s"]]
	
	nuS=Inf; bSmS=0; bS=0; locTrendFractS=0;  #these initializations are important, do not remove. 
	#t=1; irun=1
	if (SEASONALITY>1) {
		sS=rep(1,SEASONALITY+h)
	}
	
	# Initialise a matrix which contains the last level value
	yf=matrix(object$paramMeans[["lastLevel"]],nrow=NUM_OF_TRIALS, ncol=h)
	
	# For each forecasting trial
	for (irun in 1:NUM_OF_TRIALS) {
	  # Obtain the relevant parameters
		indx=sample(nrow(object$params[["l"]]),1)
		prevLevel=object$params[["lastLevel"]][indx]
		levSmS=object$params[["levSm"]][indx]
		if (!is.null(nu)) {
			nuS=nu[indx]
		} 
		
		powTrendS=object$params[["powTrend"]][indx]
		coefTrendS=object$params[["coefTrend"]][indx]
		
		if(!is.null(lastB)) {
			bS=lastB[indx]
			bSmS= object$params[["bSm"]][indx]
			locTrendFractS=object$params[["locTrendFract"]][indx]
		}
		
		sigmaS <- object$params[["sigma"]][indx]
		if (!is.null(lastSmoothedInnovSize)) {
			innovSize=lastSmoothedInnovSize[indx]
			innovSmS=object$params[["innovSm"]][indx]
			offsetsigmaS=object$params[["offsetSigma"]][indx]
		} else if (!is.null(powx)) {
			powxS=powx[indx]
			offsetsigmaS=object$params[["offsetSigma"]][indx]
		}
		
		#t=1
		if (SEASONALITY>1) {
			sS[1:SEASONALITY]=s[indx,(ncol(s)-SEASONALITY+1):ncol(s)]
			for (t in 1:h) {
				seasonA=sS[t]
				## From eq.6.a
				expVal=(prevLevel+ coefTrendS*abs(prevLevel)^powTrendS)*seasonA;
				if (!is.null(powx)) {
					omega=sigmaS*(abs(prevLevel))^powxS+offsetsigmaS
				} else if (!is.null(lastSmoothedInnovSize)) {
					omega=sigmaS*innovSize + offsetsigmaS
				} else {
					omega=sigmaS
				}
				
				# Find the t-dist error
				error=rst(n=1, xi=0 ,	omega=omega, alpha=0, nu=nuS)
				
				# Fill in the matrix of predicted value
				yf[irun,t]=min(MAX_VAL,max(MIN_VAL,expVal+error))
				
				# find the currLevel
				currLevel=max(MIN_VAL,levSmS*yf[irun,t]/seasonA + (1-levSmS)*prevLevel) ;
				
				if (currLevel>MIN_VAL) {
					prevLevel=currLevel
				} 
				sS[t+SEASONALITY] <- sS[t];
			}	#through horizons
			# yf[irun,]
		} else { #nonseasonal
			for (t in 1:h) {
				expVal=prevLevel + coefTrendS*(abs(prevLevel))^powTrendS + locTrendFractS*bS
				if (!is.null(powx)) {
					omega=sigmaS*(abs(prevLevel))^powxS+offsetsigmaS
				} else if (!is.null(lastSmoothedInnovSize)) {
					omega=sigmaS*innovSize + offsetsigmaS
				} else {
					omega=sigmaS
				}
				error=rst(n=1, xi=0, omega=omega, alpha=0, nu=nuS)
				
				yf[irun,t]=min(MAX_VAL,max(MIN_VAL, expVal+error))
				currLevel=max(MIN_VAL,levSmS*yf[irun,t] + (1-levSmS)*prevLevel) ;
				
				if (currLevel>MIN_VAL) {
					bS= bSmS*(currLevel-prevLevel)+(1-bSmS)*bS #but bSmS and bS may be==0 so then noop
					prevLevel=currLevel
				} 
			} #through horizons
			# yf[irun,]
		}
	} #through trials (simulations)
	  
  out$yf <- yf
  
  #' @importFrom stats ts
  out$mean <- ts(apply(yf, 2, mean),frequency=frequency(out$x),start=start.f   )

  #' @importFrom stats quantile
  avgYfs=apply(yf,2,quantile,probs=quantiles)
  
  out$median <- ts(avgYfs[indexOfMedian,])
	if (indexOfMedian>1) {
		out$lower <- ts(t(avgYfs[1:(indexOfMedian-1),]))
		out$upper <- ts(t(avgYfs[(indexOfMedian+1):(indexOfMedian+length(upperPercentiles)),]))	
	}
	
  out$level <- level
  
  # Assign tsp attributes for consistency
  tsp(out$median) <- tsp(out$lower) <- tsp(out$upper) <- tsp(out$mean)
	  
  class(out) <- "forecast"
  out
  
}

