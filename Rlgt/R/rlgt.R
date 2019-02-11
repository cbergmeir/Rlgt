#' @title Fit an Rlgt model
#' @description The main function to fit an rlgt model. It fits the parameter values with MCMC.
#' @param y time-series data for training (provided as a numeric vector, or a ts, or msts object).
#' @param seasonality This specification of seasonality will be overridden by frequency of y, if y is of ts or msts class. 
#' 1 by default, i.e. no seasonality.
#' @param seasonality2 Second seasonality. If larger than 1, a dual seasonality model will be used. 
#' However, this is experimental. If not specified and multiple seasonality time series (of msts class) is used,
#' a single seasonality model will be applied, one with seasonality equal to the largest of seasonalities of the time series. 
#' 1 by default, i.e. no seasonality or single seasonality.
#' @param seasonality.type Either "multiplicative" (default) or "generalized". 
#' The latter seasonality generalizes additive and multiplicative seasonality types.
#' @param error.size.method Function providing size of the error. Either "std" (monotonically, but slower than proportionally, growing with the series values) or 
#' "innov" (proportional to a smoothed abs size of innovations, i.e. surprises)  
#' @param level.method "HW",  "seasAvg", "HW_sAvg". Here, "HW" follows Holt-Winters approach. 
#' "seasAvg" calculates level as a smoothed average of the last seasonality number of points (or seasonality2 of them for the dual seasonality model),
#' and HW_sAvg is an weighted average of HW and seasAvg methods. 
#' @param xreg Optionally, a vector or matrix of external regressors, which must have the same number of rows as y. 
#' @param control list of control parameters, e.g. hyperparameter values for the model's prior distributions, number of fitting interations etc.  
#' @param verbose whether verbose information should be printed (Boolean value only), default \code{FALSE}.
#' @return \code{\link{rlgtfit}} object
#' @examples
# \dontrun{
#' # The following is a toy example that runs within a few seconds. To get good 
#' # fitting results the number of iterations should be set to at least 2000, and 
#' # 4 chains should be used (the default). To speed up computation the number of 
#' # cores should also be adjusted (default is 4).
#' 
#' rlgt_model <- rlgt(lynx, 
#'        control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=50, NUM_OF_CHAINS = 1, 
#'                             NUM_OF_CORES = 1), verbose=TRUE)
#'
#' # print the model details
#' print(rlgt_model)
#}
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
rlgt <- function(   #y=trainData; seasonality=12; seasonality2=1; seasonality.type="multiplicative"; error.size.method="std"; level.method="seasAvg";xreg = NULL; control=rlgt.control(NUM_OF_ITER=5000); verbose=TRUE; library(rstan)
	y,  
 	seasonality=1, seasonality2=1, 
 	seasonality.type=c("multiplicative","generalized"),
 	error.size.method=c("std","innov"),
 	level.method=c("HW", "seasAvg","HW_sAvg"),
 	xreg = NULL,
 	control=rlgt.control(), 
 	verbose=FALSE) {

  oldWidth=options("width")
  options(width=180)
  
  # for safety
  #model.type <- model.type[1]
  error.size.method <- error.size.method[1]
  
  seasonalityMethodId=0
  seasonality.type <- seasonality.type[1]
  if (seasonality.type=="generalized") {
	  seasonalityMethodId=1
  }
	
  levelMethodId=0
  level.method<-level.method[1]
	#yes, no levelMethodId==1, just history :-)
  if (level.method=="seasAvg") {
    levelMethodId=2
  } else if (level.method=="HW_sAvg") {
		levelMethodId=3
	}
  useSmoothingMethodForError<-error.size.method=="innov"
  nChains<-control$NUM_OF_CHAINS
  nCores<-control$NUM_OF_CORES
  addJitter<-control$ADD_JITTER
  use.regression <- !is.null(xreg)
  
  if (inherits(y,'msts')) {
    if (seasonality2<=1) {
      seasonality <- max(attributes(y)$msts)
    } else {
      seasonality <- min(attributes(y)$msts)
    }
  } else if (inherits(y,'ts')) {
    seasonality <- frequency(y)
  } 
  
  if (seasonality>1 || seasonality2>1) {
    if (seasonality>1 && seasonality2>1) { #dual seasonality
      if (seasonality>seasonality2) { #swap seasonalities
        temp <- seasonality
        seasonality <- seasonality2
        seasonality2 <- temp
      }
      model.type <- "S2GT"
    } else {#single seasonality
      if (seasonality2>1) { #swap seasonalities
        seasonality <- seasonality2
        seasonality2 <- 1
      }
      model.type <- "SGT"
    } 
  } else { #non-seasonal
    model.type <- "LGT"
  }
  
  if (seasonality <= 1 && levelMethodId != 0) {
    print("Warning: nonstandard level methods implemented only for seasonality models")
  }  
  
  model <- initModel(model.type = model.type,   #here
                     use.regression = use.regression, 
										 seasonalityMethodId=seasonalityMethodId,
										 levelMethodId=levelMethodId,  
                     useSmoothingMethodForError=useSmoothingMethodForError)
  
  MAX_RHAT_ALLOWED <- control$MAX_RHAT_ALLOWED
  NUM_OF_ITER <- control$NUM_OF_ITER
  MAX_NUM_OF_REPEATS <- control$MAX_NUM_OF_REPEATS
  
  if (nCores>1) {
    rstan_options(auto_write = TRUE)
    options(mc.cores = nCores)    
  }
  
  y_org <- y
  y <- as.numeric(y)
  n <- length(y)
  
  #' @importFrom stats rnorm
  if(addJitter) y <- y + rnorm(n, 0, sd = abs(min(y)) * 0.0001)
  
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
							 SEASONALITY_TYPE=seasonalityMethodId,
               USE_REGRESSION=as.integer(use.regression), 
               LEVEL_CALC_METHOD=levelMethodId,
               #if this is a regression call, these three values will be overwritten in a moment below, 
               # I can't make it J==1, becasue then REG_CAUCHY_SD becomes number, not vector
               J=2, 
               xreg=matrix(0,nrow=n, ncol=2), 
               REG_CAUCHY_SD=rep(1,2)) 
  
  if (use.regression) {
    if (!is.matrix(xreg)) { # convert non-matrix to matrix
      xreg <- as.matrix(xreg)
    }
    if (nrow(xreg) != n) {
      stop("Error: Number of rows(", nrow(xreg),") supplied in xreg != length of y(", n, ").")
    }
    data[['xreg']] <- xreg
    data[['J']]    <- ncol(xreg)
    regCauchySd <- mean(y)/apply(xreg,2,mean)/control$CAUCHY_SD_DIV #vector
    if (ncol(xreg)==1) dim(regCauchySd)=1
    data[['REG_CAUCHY_SD']] <- regCauchySd
  } 
  
	avgDiff=mean(abs(diff(y)))
	
	
  ### Repeat until Rhat is in an acceptable range (i.e. convergence is reached)
  avgRHat <- 1e200
  irep <- 1
  for (irep in 1:MAX_NUM_OF_REPEATS) {
    initializations <- list()
    for (irr in 1:nChains) {
      initializations[[irr]] <- list( 
        innovSizeInit = abs(rnorm(1, 0, CauchySd)),#used only for *GTe models 
        bInit=rnorm(1,mean=0, sd=y[1]/control$CAUCHY_SD_DIV), 
        coefTrend=rnorm(1,mean=0, sd=y[1]/control$CAUCHY_SD_DIV),
        sigma=runif(1,min=avgDiff/100, max=avgDiff/10),
        offsetSigma=runif(1,min=avgDiff/100, max=avgDiff/10)
      )
      if (use.regression) {
        initializations[[irr]][['regCoef']] <- rnorm(ncol(xreg),mean=0, sd=regCauchySd)
        if (ncol(xreg)==1) dim(initializations[[irr]][['regCoef']])=1
        initializations[[irr]][['regOffset']] <- rnorm(1,mean=0, sd=mean(regCauchySd))
      }
      if (seasonalityMethodId==0 || seasonalityMethodId==2) {
				initializations[[irr]][['initS']] <- rnorm(seasonality, 0, 0.05) #exp()
				initializations[[irr]][['initS2']] <- rnorm(seasonality2, 0, 0.05)
      } else {  #generalized
				initializations[[irr]][['initS']] <- rnorm(seasonality, mean=0, sd=min(y[1:seasonality])*0.003) 
				initializations[[irr]][['initS2']] <- rnorm(seasonality2, mean=0, sd=min(y[1:seasonality2])*0.003)
      }
    }
    
    #Double the number of iterations for the next cycle
    numOfIters <- NUM_OF_ITER * 2 ^ (irep - 1)
    samples = rstan::sampling(
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
    ainfo <- summary(samples)
    RHats <- ainfo$summary[,10]
    RHats <- as.numeric(RHats[is.finite(RHats)])
    currRHat <- mean(RHats, na.rm = T)
    if (currRHat <= MAX_RHAT_ALLOWED) {
      avgRHat <- currRHat
      if(verbose) print(samples)
      print(paste("avgRHat",avgRHat))
      break
    } else {
      if (currRHat<avgRHat) {#in the worst case this is at least once executed, because avgRHat is initialized high
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
    print(summary(do.call(rbind, args = get_sampler_params(samples, inc_warmup = F)), digits = 2)) #diagnostics including step sizes and tree depths
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
	if (levelMethodId==3) {
		params[["lastLevel0"]] <- params[["l0"]][,ncol(params[["l0"]])]
	}
	
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
			seasonalityMethodId=seasonalityMethodId, levelMethodId=levelMethodId,  
      useSmoothingMethodForError=useSmoothingMethodForError,
      seasonality=seasonality, seasonality2=seasonality2,   
      model, params, control, samples)
  out
}
