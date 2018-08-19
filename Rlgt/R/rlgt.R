#' @title Fit an Rlgt model
#' @description The main function to fit an rlgt model. It fitted the parameter values through MCMC simulations.
#' @param y time-series data for training (provided as a vector or a ts object).
#' @param model.type a chosen model from the Rlgt package.
#' @param xreg Optionally, a vector or matrix of external regressors, which must have the same number of rows as y.
#' @param control list of control parameters, i.e. hyperparameter values for the model's prior distribution.  
#' @param nChains number of MCMC chains ro be produced. The number must be >=1. 
#' More chains will improve the quality of the model, but will significantly increase 
#' the run-time of the MCMC sampling procedure.Perhaps, an optimal number is around 4.
#' @param nCores number of cores to be used. For performance reasons, it should be set to equal 
#' to nChains, in case of nChains <= the number of cores available on the computer.  
#' @param addJitter adding a bit of jitter is helping Stan in case of some flat time-series data.
#' @param verbose whether verbose information should be printed (true/false value only)
#' @return rlgtfit object
#' @examples
#' \dontrun{
#' rlgt_model <- rlgt(lynx, model="LGT", nCores=4, nChains=4,
#' control=lgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=2000), 
#' verbose=TRUE)

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
rlgt <- function(y, model.type=c("LGT", "LGTe", "SGT", "SGTe", "S2GT", "gSGT", "Trend"), 
                 xreg = NULL,
                 control=rlgt.control(), nChains=2, nCores=2, 
                 addJitter=TRUE, verbose=FALSE) {
  # for safety
  model.type <- model.type[1]
  modelIsSeasonal <- model.type %in% c("SGT", "S2GT", "SGTe", "gSGT")
  has.regression <- !is.null(xreg)
  use.regression <- FALSE
  if (has.regression) {
    # regression components detected
    if (model.type %in% c("LGT", "SGT")) {
      use.regression <- TRUE
    } else {
      message("Current model do not support regression. Regression variables will be ignored.")
    }
  }

  if(!inherits(model.type, "RlgtStanModel")) {
    model <- initModel(model.type = model.type,
                       use.regression = use.regression)
  }
  
  MAX_RHAT_ALLOWED <- control$MAX_RHAT_ALLOWED
  NUM_OF_ITER <- control$NUM_OF_ITER
  MAX_NUM_OF_REPEATS <- control$MAX_NUM_OF_REPEATS
  
  if (nCores>1) {
    rstan_options(auto_write = TRUE)
    options(mc.cores = nCores)    
  }
  
  y.orig <- y
  n <- length(y)
  
  #' @importFrom stats rnorm
  if(addJitter) y <- y + rnorm(n, 0, sd=abs(min(y)) * 0.0001)
  
  CauchySd <- max(y) / control$CAUCHY_SD_DIV
  
  SEASONALITY <- control$SEASONALITY
  SEASONALITY2 <- control$SEASONALITY2
  
  #' @importFrom stats frequency
  if (SEASONALITY <= 1 && frequency(y) > 1 && modelIsSeasonal) {
    SEASONALITY <- frequency(y)
    print(paste0("Seasonality not specified, but the data is seasonal. Inferring seasonality equal to ",SEASONALITY))
    control$SEASONALITY <- SEASONALITY  #will be used in forecast()
  }
  
  data <- list(CAUCHY_SD=CauchySd, 
               SEASONALITY=SEASONALITY, 
               SEASONALITY2=SEASONALITY2,
               MIN_POW_TREND=control$MIN_POW_TREND, 
               MAX_POW_TREND=control$MAX_POW_TREND, 
               MIN_SIGMA=control$MIN_SIGMA,  
               MIN_NU=control$MIN_NU,  
               MAX_NU=control$MAX_NU, 
               POW_SEASON_ALPHA=control$POW_SEASON_ALPHA, 
               POW_SEASON_BETA=control$POW_SEASON_BETA,
               POW_TREND_ALPHA=control$POW_TREND_ALPHA, 
               POW_TREND_BETA=control$POW_TREND_BETA,
               y=y, N=n, SKEW=control$SKEW) # to be passed on to Stan
  
  if(use.regression){
    data[['xreg']] <- xreg
    data[['J']] <- ncol(xreg)
    if(nrow(xreg) != n){
      stop("Error: Number of rows supplied in regression matrix does not match length of y!")
    }
  }
  
  ### Repeat until Rhat is in an acceptable range (i.e. convergence is reached)
  avgRHat <- 1e200
  irep <- 1
  for (irep in 1:MAX_NUM_OF_REPEATS) {
    initializations <- list()
    for (irr in 1:nChains) {
      initializations[[irr]]=list( 
        ### Initialise seasonality factors
        # for non-seasonal models it is not necessary, 
        # but makes code simpler and is not a big overhead
        initSu = rnorm(SEASONALITY, 1, 0.05), 
        initSu2 = rnorm(SEASONALITY2, 1, 0.05),
        innovSizeInit = abs(rnorm(1, 0, CauchySd))#used only for *GTe models 
      )
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
    
    # Find the mean, but for ETS components average based on each t
    #paramMeans[[param]] <- if(param %in% c("l", "b", "s")) apply(params[[param]],2,mean) else mean(params[[param]])
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
  
  out <- rlgtfit(y.orig, model.type, use.regression = use.regression,
             model, params, control, samples)
  out
}