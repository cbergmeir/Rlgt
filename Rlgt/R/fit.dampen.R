
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
#' @returnType dampen
#' @return dampenmodel
#' @examples
#'\dontrun{
#' lgt_model <- fit.lgt(lynx, model="LGT", nCores=4, nChains=4,
#' control=lgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=2000), 
#' verbose=TRUE)

#' # print the model details
#' print(lgt_model)
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
#' 
fit.dampen <- function(y, model= c("Dampen", "SDampen", "TDampen", "TSDampen"), 
  control=lgt.control(), nChains=2, nCores=2, addJitter=TRUE, verbose=FALSE) {
  
  if(!inherits(model, "RlgtStanModelDampen")) {
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
    MIN_SIGMA=control$MIN_SIGMA, MIN_NU=control$MIN_NU,  
    MAX_NU=control$MAX_NU,
    POW_SIGMA_ALPHA=control$POW_SIGMA_ALPHA, 
    POW_SIGMA_BETA=control$POW_SIGMA_BETA,
    y=y, N=n, SKEW=control$SKEW) # to be passed on to Stan
  
  ### Repeat until Rhat is in an acceptable range (i.e. convergence is reached)
  avgRHat=1e200; irep=1
  for (irep in 1:MAX_NUM_OF_REPEATS) {
    print("irep:", irep)
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

  
  
  out <- dampen(y.orig, model, params, paramMeans, SEASONALITY, samples)
  
  out
  
}