#' @title Fit an Rlgt model
#' @description The main function to fit an rlgt model. It fits the parameter values with MCMC.
#' @param y time-series data for training (provided as a numeric vector, or a ts, or msts object).
#' @param seasonality This specification of seasonality will be overridden by frequency of y, if y is of ts or msts class. 
#' 1 by default, i.e. no seasonality.
#' @param seasonality2 Second seasonality. If larger than 1, a dual seasonality model will be used. 
#' This specification of seasonality will be overridden by the second seasonality of y, if y is of msts class. 
#' 1 by default, i.e. no seasonality or simgle seasonality.
#' @param seasonality.type Either "multiplicative" (default) or "generalized". 
#' The latter seasonality generalizes additive and multiplicative seasonality types.
#' @param error.size.method It chooses a function providing size of the error. Either "std" (monotonically, but slower than proportionally, growing with the series values) or 
#' "innov" (proportional to a smoothed abs size of innovations, i.e. surprises)  
#' @param level.method one of "classical", "seasAvg", or "seas2Avg". Here, "classical" is normal Holt-Winters, with level divided by seasonality. "seasAvg" is with a window with the size of the seasonality, we take the median of the full seasonality, i.e., not y divided by the seasonality. "seas2Avg" is the same as seasAvg, taking as window the size of the larger of the two seasonalities.
#' @param xreg Optionally, a vector or matrix of external regressors, which must have the same number of rows as y. Default is "std".
#' @param control list of control parameters, e.g. hyperparameter values for the model's prior distributions, number of fitting interations etc.  
#' @param verbose whether verbose information should be printed (Boolean value only), default \code{FALSE}.
#' @return \code{\link{rlgtfit}} object
#' @examples
#' \dontrun{
#' rlgt_model <- rlgt(lynx, control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=2000), 
#'      verbose=TRUE)
#'
#' # print the model details
#' print(rlgt_model)
#'}
#'
#'\dontrun{demo(exampleScript)}
#'
#' @importFrom rstan rstan_options
#' @importFrom rstan sampling
#' @importFrom rstan get_sampler_params
#' @importFrom rstan extract
#' @importMethodsFrom rstan summary
#' @importFrom sn rst
#' @export
rlgt <- function(y,  
						seasonality=1, seasonality2=1, 
						seasonality.type=c("multiplicative","generalized"),
						error.size.method=c("std","innov"),
						level.method=c("classical","seasAvg"),
            xreg = NULL,
            control=rlgt.control(), 
						verbose=FALSE) {
					
	oldWidth=options("width")
	options(width=180)
	
  # for safety
  #model.type <- model.type[1]
	error.size.method <- error.size.method[1]
	seasonality.type <- seasonality.type[1]
	level.method<-level.method[1]
	levelMethodId=0
	if (level.method=="seasAvg") {
		levelMethodId=1
	} 
	useGeneralizedSeasonality<-seasonality.type=="generalized"
	useSmoothingMethodForError<-error.size.method=="innov"
	nChains<-control$NUM_OF_CHAINS
	nCores<-control$NUM_OF_CORES
	addJitter<-control$ADD_JITTER
	use.regression <- !is.null(xreg)
	
	if (inherits(y,'msts')) {
		seasonality=attributes(y)$msts[1]
		seasonality2=attributes(y)$msts[2]
	} else if (inherits(y,'ts')) {
		seasonality=frequency(y)
	} 
	
	if (seasonality>1 || seasonality2>1) {
		if (seasonality>1 && seasonality2>1) { #dual seasonality
			if (seasonality>seasonality2) { #swap seasonalities
				temp=seasonality
				seasonality=seasonality2
				seasonality2=temp
			}
			model.type="S2GT"
		} else {#single seasonality
			if (seasonality2>1) { #swap seasonalities
				seasonality=seasonality2
				seasonality2=1
			}
			model.type="SGT"
		} 
	} else { #non-seasonal
			model.type="LGT"
	}
		
	if (seasonality2<=1 && levelMethodId!=0) {
		print("Warning: nonstandard level methods implemented only for dual seasonality models. level.method will be ignored")
	}  
	
  model <- initModel(model.type = model.type,   #here
              use.regression = use.regression, 
							useGeneralizedSeasonality=useGeneralizedSeasonality,
							useSmoothingMethodForError=useSmoothingMethodForError)
  
  MAX_RHAT_ALLOWED <- control$MAX_RHAT_ALLOWED
  NUM_OF_ITER <- control$NUM_OF_ITER
  MAX_NUM_OF_REPEATS <- control$MAX_NUM_OF_REPEATS
  
  if (nCores>1) {
    rstan_options(auto_write = TRUE)
    options(mc.cores = nCores)    
  }
	
  y_org <- y
	y <-as.numeric(y)
  n <- length(y)
  
  #' @importFrom stats rnorm
  if(addJitter) y <- y + rnorm(n, 0, sd=abs(min(y)) * 0.0001)
  
  CauchySd <- max(y) / control$CAUCHY_SD_DIV
	
	
  
  data <- list(CAUCHY_SD=CauchySd, 
               SEASONALITY=as.integer(seasonality),
							 SEASONALITY_F=seasonality,
               SEASONALITY2=as.integer(seasonality2),
							 SEASONALITY2_F=seasonality2,
							 NUM_OF_SEASON_INIT_CYCLES=control$NUM_OF_SEASON_INIT_CYCLES,
               MIN_POW_TREND=control$MIN_POW_TREND, 
               MAX_POW_TREND=control$MAX_POW_TREND, 
               MIN_SIGMA=control$MIN_SIGMA,  
               MIN_NU=control$MIN_NU,  
               MAX_NU=control$MAX_NU, 
               POW_SEASON_ALPHA=control$POW_SEASON_ALPHA, 
               POW_SEASON_BETA=control$POW_SEASON_BETA,
               POW_TREND_ALPHA=control$POW_TREND_ALPHA, 
               POW_TREND_BETA=control$POW_TREND_BETA,
               y=y, N=n, 
							 USE_SMOOTHED_ERROR=as.integer(useSmoothingMethodForError),
					 		 USE_GENERALIZED_SEASONALITY=as.integer(useGeneralizedSeasonality),
							 USE_REGRESSION=as.integer(use.regression), 
							 LEVEL_CALC_METHOD=levelMethodId,
							 J=2, xreg=matrix(0,nrow=n, ncol=2), REG_CAUCHY_SD=rep(1,2)) #if this is a regression call, these three values will be overwritten in a moment below, I can't make it J==1, becasue then REG_CAUCHY_SD becomes number, not vector
  
  if (use.regression) {
    if (is.null(dim(xreg))) { # convert non-matrix to matrix
      xreg <- as.matrix(xreg)
    }
    if (nrow(xreg) != n) {
      stop("Error: Number of rows supplied in regression matrix does not match length of y!")
    }
    data[['xreg']] <- xreg
    data[['J']]    <- ncol(xreg)
		regCauchySd <- mean(y)/apply(xreg,2,mean)/control$CAUCHY_SD_DIV #vector
		data[['REG_CAUCHY_SD']] <- regCauchySd
  } 
	
	
  
  ### Repeat until Rhat is in an acceptable range (i.e. convergence is reached)
  avgRHat <- 1e200
  irep <- 1
  for (irep in 1:MAX_NUM_OF_REPEATS) {
    initializations <- list()
    for (irr in 1:nChains) {
      initializations[[irr]]=list( 
        # These initializations are strictly not necessary, but they improve performance (accuracy and calculation speed)
        innovSizeInit = abs(rnorm(1, 0, CauchySd)),#used only for *GTe models 
				bInit=rnorm(1,mean=0, sd=y[1]/control$CAUCHY_SD_DIV),
				coefTrend=rnorm(1,mean=0, sd=0.01),
				sigma=runif(1,min=control$MIN_SIGMA, max=0.01),
				offsetSigma=runif(1,min=control$MIN_SIGMA, max=0.01)
      )
			if (use.regression) {
				initializations[[irr]][['regCoef']] <- rnorm(ncol(xreg),mean=0, sd=regCauchySd)
				dim(initializations[[irr]][['regCoef']]) <- ncol(xreg)
				initializations[[irr]][['regOffset']] <- rnorm(1,mean=0, sd=mean(regCauchySd))
				dim(initializations[[irr]][['regOffset']]) <- 1
			}
			if (useGeneralizedSeasonality) {
				initializations[[irr]][['initSu']] <- rnorm(seasonality, mean=0, sd=min(y[1:seasonality])*0.003) 
				initializations[[irr]][['initSu2']] <- rnorm(seasonality2, mean=0, sd=min(y[1:seasonality2])*0.003)
			} else {  #multiplicative
				initializations[[irr]][['initSu']] <- rnorm(seasonality, 1, 0.05) 
				initializations[[irr]][['initSu2']] <- rnorm(seasonality2, 1, 0.05)
			}
    }
    
    #Double the number of iterations for the next cycle
    numOfIters <- NUM_OF_ITER * 2 ^ (irep - 1)
    samples1 <-
      rstan::sampling(
        control = list(adapt_delta = control$ADAPT_DELTA, 
                       max_treedepth = control$MAX_TREE_DEPTH),
        object = model$model,   
        data = data, 
        init = initializations,
        pars = model$parameters,
        iter = numOfIters,
        chains = nChains,
        cores = nCores,
        open_progress = F,
        refresh = if(verbose) numOfIters / 5 else - 1)
    
    ### Get the Rhat values
    ainfo <- summary(samples1)
    RHats <- ainfo$summary[,10]
    RHats <- as.numeric(RHats[is.finite(RHats)])
    currRHat <- mean(RHats, na.rm = T)
    if (currRHat <= MAX_RHAT_ALLOWED) {
      samples <- samples1
      avgRHat <- currRHat
      if(verbose) print(samples)
      print(paste("avgRHat",avgRHat))
      break
    } else {
      if (currRHat<avgRHat) {#in the worst case this is at least once executed, because avgRHat is initialized high
        samples <- samples1
        avgRHat <- currRHat
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
  
  if(verbose) {
    print(summary(do.call(rbind, 
                          args = get_sampler_params(samples1, 
                                                    inc_warmup = F)), 
                  digits = 2)) #diagnostics including step sizes and tree depths
  }
  
  params <- list()
  #paramMeans <- list()
  
  # Extract all of the parameter means
  for(param in model$parameters) {
    params[[param]] <- extract(samples)[[param]]
  }
  
  #special processing for vector params, where only the last value counts
  # This includes level, trend, and innov size
  params[["lastLevel"]] <- params[["l"]][,ncol(params[["l"]])]
  #paramMeans[["lastLevel"]] <- mean(params[["lastLevel"]])
  
  if ("b" %in% model$parameters) {
    params[["lastB"]] <- params[["b"]][,ncol(params[["b"]])]
    #paramMeans[["lastB"]] <- mean(params[["lastB"]])
  }
  if ("smoothedInnovSize" %in% model$parameters) {
    params[["lastSmoothedInnovSize"]] <- params[["smoothedInnovSize"]][,ncol(params[["smoothedInnovSize"]])]
    #paramMeans[["lastSmoothedInnovSize"]] <- mean(params[["lastSmoothedInnovSize"]])
  }
  
	options(width=oldWidth[[1]])
	#correct: y_org. Fitting has already been done. y_org is either numeric, or ts, or msts 
  out <- rlgtfit(y_org, model.type, use.regression = use.regression, 
			useGeneralizedSeasonality=useGeneralizedSeasonality, levelMethodId=levelMethodId,  
			useSmoothingMethodForError=useSmoothingMethodForError,
			seasonality=seasonality, seasonality2=seasonality2,   
      model, params, control, samples)
  out
}
