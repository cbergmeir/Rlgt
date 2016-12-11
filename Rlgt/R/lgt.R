# Author: Slawek Smyl 
# Oct 2016
# An example fitting the (non-seasonal) LGT to one of the M3 yearly curves
###############################################################################
# install.packages("Mcomp")
# install.packages("rstan") #but first install RTools


#' Initialize a stan model that uses the (non-seasonal) LGT
#' 
#' @title Initialize a non-seasonal LGT stan model
#' @param silent output stan compiler messages yes/no
#' @returnType 
#' @return 
#' 
#' @importFrom rstan stan_model
#' @export
init.lgt <- function(silent=TRUE){
  
#non-seasonal model
  modString = '
      data {  
      real<lower=0> CAUCHY_SD;
      real MIN_POW;  real MAX_POW;
      real<lower=0> MIN_SIGMA;
      real<lower=1> MIN_NU; real<lower=1> MAX_NU;
      int<lower=1> N;
      vector<lower=0>[N] y;
      }
      parameters {
      real<lower=MIN_NU,upper=MAX_NU> nu; 
      real<lower=0> sigma;
      real <lower=0,upper=1>levSm;
      real <lower=0,upper=1> bSm;
      real <lower=0,upper=1> powx;
      real <lower=0,upper=1> powTrendBeta;
      real coefTrend;
      real <lower=MIN_SIGMA> offsetSigma;
      real <lower=-1,upper=1> locTrendFract;
      } 
      transformed parameters {
      real <lower=MIN_POW,upper=MAX_POW>powTrend;
      vector[N] l; vector[N] b;
      
      l[1] = y[1]; b[1] = 0;
      powTrend= (MAX_POW-MIN_POW)*powTrendBeta+MIN_POW;
      
      for (t in 2:N) {
      l[t]  = levSm*y[t] + (1-levSm)*l[t-1] ;
      b[t]  = bSm*(l[t]-l[t-1]) + (1-bSm)*b[t-1] ;
      }
      }
      model {
      sigma ~ cauchy(0,CAUCHY_SD) T[0,];
      offsetSigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];
      coefTrend ~ cauchy(0,CAUCHY_SD);
      
      for (t in 2:N) {
      y[t] ~ student_t(nu, l[t-1]+coefTrend*fabs(l[t-1])^powTrend+locTrendFract*b[t-1], 
      sigma*fabs(l[t-1])^powx+ offsetSigma);
      }
      }
      '  
  
  stanModel <- list()
  
  stanModel[["parameters"]] <- c("l", "b", "nu", "sigma", "levSm",  "bSm", 
      "powx", "coefTrend",  "powTrend", "offsetSigma", "locTrendFract")
  
  out <- capture.output(stanModel[["model"]] <- stan_model(model_code=modString))
  #stanModel[["mod_output"]]
  
  if(!silent) print(out)
  
  class(stanModel) <- c("stanModel", "stanModelLGT")
  stanModel  
  
} 


#' Initialize a stan model that is a trend-only, Gaussian error, homoscedastic version of LGT
#' 
#' @title Initialize a stan model that uses trend only
#' @param silent output stan compiler messages yes/no
#' @returnType 
#' @return 
#' 
#' @importFrom rstan stan_model
#' @export
init.trend <- function(silent=TRUE){
  
  modString = '
      data {  
      real<lower=0> CAUCHY_SD;
      real MIN_POW;  real MAX_POW;
      real<lower=0> MIN_SIGMA;
      //real<lower=1> MIN_NU; real<lower=1> MAX_NU;
      int<lower=1> N;
      vector<lower=0>[N] y;
      }
      parameters {
      real<lower=MIN_SIGMA> sigma;
      real <lower=0,upper=1>levSm;
      real <lower=0,upper=1> bSm;
      real <lower=0,upper=1> powTrendBeta;
      real coefTrend;
      real <lower=-1,upper=1> locTrendFract;
      } 
      transformed parameters {
      real <lower=MIN_POW,upper=MAX_POW>powTrend;
      vector[N] l; vector[N] b;
      
      l[1] = y[1]; b[1] = 0;
      powTrend= (MAX_POW-MIN_POW)*powTrendBeta+MIN_POW;
      
      for (t in 2:N) {
      l[t]  = levSm*y[t] + (1-levSm)*l[t-1] ;
      b[t]  = bSm*(l[t]-l[t-1]) + (1-bSm)*b[t-1] ;
      }
      }
      model {
      sigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];
      coefTrend ~ cauchy(0,CAUCHY_SD);
      
      for (t in 2:N) {
      y[t] ~ normal(l[t-1]+coefTrend*fabs(l[t-1])^powTrend+locTrendFract*b[t-1], sigma);
      }
      }
      '   
  
  stanModel <- list()
  
  stanModel[["parameters"]] <- c("l", "b", "sigma", "levSm",  "bSm", "coefTrend", 
      "powTrend", "locTrendFract")
  
  out <- capture.output(stanModel[["model"]] <- stan_model(model_code=modString))
  #stanModel[["mod_output"]]
  
  if(!silent) print(out)
  
  class(stanModel) <- c("stanModel", "stanModelTrend")
  stanModel  
  
} 









#myMod <- initPredictionModel.lgt()
#length(myMod$output)



#' This is a function that initializes and sets the parameters of the algorithm. 
#' It generates a list of parameters, to be used with the \code{\link{fit.lgt}} function. 
#'
#' @title Sets and initializes the main parameters of the algorithm
#' @param MAX_RHAT_ALLOWED 
#' @param NUM_OF_ITER 
#' @param MAX_NUM_OF_REPEATS 
#' @param CAUCHY_SD_DIV 
#' @param MIN_SIGMA 
#' @param MIN_NU 
#' @param MAX_NU 
#' @param MIN_POW 
#' @param MAX_POW 
#' @returnType 
#' @return 
#' @export
lgt.control <- function(
    
    MAX_RHAT_ALLOWED=1.005, #see Stan's manual
    NUM_OF_ITER=5000,
    MAX_NUM_OF_REPEATS=3,
    CAUCHY_SD_DIV=200,
    MIN_SIGMA=0.001, # for numerical stability
    MIN_NU=1.5,
    MAX_NU=20,
    MIN_POW=-0.5,
    MAX_POW=1) {
  
  list(CAUCHY_SD_DIV=CAUCHY_SD_DIV,
      MAX_RHAT_ALLOWED=MAX_RHAT_ALLOWED,
      NUM_OF_ITER=NUM_OF_ITER,
      MAX_NUM_OF_REPEATS=MAX_NUM_OF_REPEATS,
      MIN_SIGMA=MIN_SIGMA,
      MIN_NU=MIN_NU,
      MAX_NU=MAX_NU,
      MIN_POW=MIN_POW,
      MAX_POW=MAX_POW)
}






#' Runs the model fitting
#' 
#' @title Runs the model fitting
#' @param y the time series
#' @param h prediction horizon
#' @param stanModel a stan model
#' @param control control arguments list
#' @param ncores number of cores
#' @param addJitter adding a bit of jitter is helping Stan in case of some flat series
#' @param verbose print verbose information yes/no
#' @returnType 
#' @return 
#' 
#' @importFrom rstan rstan_options
#' @importFrom rstan sampling
#' @importFrom rstan get_sampler_params
#' @importFrom rstan extract
#' @importMethodsFrom rstan summary
#' @importFrom sn rst
#' @export
fit.lgt <- function(y, stanModel=NULL, control=lgt.control(), ncores=1, addJitter=TRUE, verbose=FALSE) {
  
  if(is.null(stanModel)) {
    warning("No model was provided, generating an LGT model. This is slow, so if you run this function on many series better generate the model once and then supply it as a parameter")
    stanModel <- init.lgt()
  }
  
  MAX_RHAT_ALLOWED=control$MAX_RHAT_ALLOWED
  NUM_OF_ITER=control$NUM_OF_ITER
  MAX_NUM_OF_REPEATS=control$MAX_NUM_OF_REPEATS
  
  STAN_CLUSTER_SIZE=ncores
  
  if(ncores>1) {
    rstan_options(auto_write = TRUE)
    options(mc.cores = STAN_CLUSTER_SIZE)    
  }
  
  y.orig <- y
  n=length(y)
  
  if(addJitter) y <- y + rnorm(n,0,sd=abs(min(y))*0.0001)
  
  
  CauchySd=max(y)/control$CAUCHY_SD_DIV
  data = list(CAUCHY_SD=CauchySd,
      MIN_SIGMA=control$MIN_SIGMA,
      MIN_NU=control$MIN_NU,
      MAX_NU=control$MAX_NU,
      MIN_POW=control$MIN_POW,
      MAX_POW=control$MAX_POW,
      y=y, N=n) # to be passed on to Stan  
  
  avgRHat=1e20; irep=1
  for (irep in 1:MAX_NUM_OF_REPEATS) {
    #initializations = list();

    samples1=
        sampling(#control=list(adapt_delta = 0.9, max_treedepth=11),
            stanModel$model,   
            data=data, 
            #init=initializations,
            pars=stanModel$parameters,
            iter=NUM_OF_ITER*2^(irep-1),
            chains=STAN_CLUSTER_SIZE,
            cores=STAN_CLUSTER_SIZE,
            open_progress=F,
            refresh = if(verbose) 1000 else -1)

    #samples1
    
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
      print (paste("trying to do better..."))
    }
    #str(samples, max.level =4)
  }#repeat if needed
  
  if(verbose) print(summary(do.call(rbind, args = get_sampler_params(samples1, inc_warmup = F)), digits = 2)) #diagnostics
  
#  stanModel[["parameters"]] <- c("l", "b", "nu", "sigma", "levSm",  "bSm", 
#      "powx", "coefTrend",  "powTrend", "offsetSigma", "locTrendFract")
  
  params <- list()
  paramMeans <- list()
  
  for(param in stanModel$parameters) {
    params[[param]] <- extract(samples)[[param]]
    paramMeans[[param]] <- if(param %in% c("l", "b")) apply(params[[param]],2,mean) else mean(params[[param]])
  }
  
  params[["lastLevel"]] <- params[["l"]][,ncol(params[["l"]])]
  params[["lastB"]] <- params[["b"]][,ncol(params[["b"]])]
  
  paramMeans[["lastLevel"]] <- mean(params[["lastLevel"]])
  paramMeans[["lastB"]] <- mean(params[["lastB"]])
  
  out <- list()

  out[["x"]] <- y.orig
  out[["stanModel"]] <- stanModel
  out[["params"]] <- params
  out[["paramMeans"]] <- paramMeans
  class(out) <- "lgt"
  
  out

}





#' This function produces forecasts from a model
#' 
#' @title produce forecasts
#' @param object 
#' @param h 
#' @param NUM_OF_TRIALS 
#' @param MIN_VAL 
#' @param MAX_VAL 
#' @param ... 
#' @returnType 
#' @return returns a forecast object compatible with the forecast package
#' @S3method forecast lgt
# @method forecast lgt
#' @importFrom forecast forecast 
#' @author bergmeir
#' @export
forecast.lgt <- function(object, h=ifelse(frequency(object$x)>1, 2*frequency(object$x), 10), NUM_OF_TRIALS=5000, 
    MIN_VAL=0.001, MAX_VAL=1e38, ...) {

  #object <- mod[["lgt"]]
  
  out <- list(model=object,x=object$x)
  #out <- object
  tspx <- tsp(out$x)
  
  if(!is.null(tspx))
    start.f <- tspx[2] + 1/frequency(out$x)
  else
    start.f <- length(out$x)+1
  
  isLGT <- inherits(object$stanModel, "stanModelLGT")
  
#  t=1; irun=1
  yf=matrix(object$paramMeans[["lastLevel"]],nrow=NUM_OF_TRIALS, ncol=h)
  
  for (irun in 1:NUM_OF_TRIALS) {
    
    indx=sample(nrow(object$params[["l"]]),1)

#TODO: I don't understand what is done here...calculating a mean over single numbers?
#    prevLevel=mean(lastLevel[indx], na.rm=T)
#    powTrendM=mean(powTrend[indx], na.rm=T)
#    coefTrendM=mean(coefTrend[indx], na.rm=T)
#    bM=mean(lastB[indx], na.rm=T)
#    levSmM=mean(levSm[indx], na.rm=T)
#    bSmM=mean(bSm[indx], na.rm=T)
#    sigmaM=mean(sigma[indx], na.rm=T)
#    locTrendFractM=mean(locTrendFract[indx], na.rm=T)
    
    prevLevel <- object$params[["lastLevel"]][indx]
    powTrendM <- object$params[["powTrend"]][indx]
    coefTrendM <- object$params[["coefTrend"]][indx]
    bM <- object$params[["lastB"]][indx]
    levSmM <- object$params[["levSm"]][indx]
    bSmM <- object$params[["bSm"]][indx]
    sigmaM <- object$params[["sigma"]][indx]
    locTrendFractM <- object$params[["locTrendFract"]][indx]
       
    
    if(isLGT) {
#      powxM=mean(powx[indx], na.rm=T)
#      nuM=mean(nu[indx], na.rm=T)
#      offsetSigmaM=mean(offsetSigma[indx], na.rm=T)
      
      powxM <- object$params[["powx"]][indx]
      nuM <- object$params[["nu"]][indx]
      offsetSigmaM <- object$params[["offsetSigma"]][indx]
    }
    
    for (t in 1:h) { 
      
      if(isLGT) {
        error=rst(n=1, 
            xi=coefTrendM*(abs(prevLevel))^powTrendM + locTrendFractM*bM , 
            omega=sigmaM*(abs(prevLevel))^powxM+offsetSigmaM, alpha=0, nu=nuM)
        yf[irun,t]=min(MAX_VAL,max(MIN_VAL,prevLevel+error))
        
      } else {
        expVal=prevLevel+coefTrendM*(abs(prevLevel))^powTrendM + locTrendFractM*bM
        error=rnorm(n=1, mean=0, sd=sigmaM)
        yf[irun,t]=min(MAX_VAL,max(MIN_VAL,expVal+error))        
      }     
      
      currLevel=max(MIN_VAL,levSmM*yf[irun,t] + (1-levSmM)*prevLevel) ;
      
      if (currLevel>MIN_VAL) {
        bM= bSmM*(currLevel-prevLevel)+(1-bSmM)*bM
        prevLevel=currLevel
      } 
    }
    # yf[irun,]
  } #through trials
  
  out$yf <- yf
  
  
  out$mean <- ts(apply(yf, 2, mean),frequency=frequency(out$x),start=start.f)

  quants=c(0.05, 0.2, 0.5, 0.8, 0.95)
  avgYfs=apply(yf,2,quantile,probs=quants)
  
  out$median <- ts(avgYfs[3,])
  out$lower <- ts(t(avgYfs[c(1,2),]))
  out$upper <- ts(t(avgYfs[c(4,5),]))
  out$level <- c(80,95)
  
#  out$median <- ts(apply(object$yf, 2, median))
#  out$lower <- ts(apply(object$yf, 2, min))
#  out$upper <- ts(apply(object$yf, 2, max))
#  out$level <- 100
  
  tsp(out$median) <- tsp(out$lower) <- tsp(out$upper) <- tsp(out$mean)
  
  class(out) <- "forecast"
  out
  
}


#' Print out some characteristics of a \code{\link{lgt}} model.
#'  
#' @title Generic print function for lgt models
#' @param x the \code{\link{lgt}} model
#' @param ... additional function parameters (currently not used)
#' @export
#' @S3method print lgt
# @method print lgt
# @rdname lgt
print.lgt <- function(x, ...) {
  if(!inherits(x, "lgt")) stop("not a legitimate lgt result")
  
 
#    lastDetails=paste("coefT=",round(coefTrendM,2), ", powT=",round(powTrendM,2),
#      ", sigma=",round(sigmaM,2),", powx=",round(powxM,2),", offsetS=",round(offsetSigmaM,2), 
#      ", nu=",round(nuM,2), ", lSm=",round(levSmM,2), 
#      ", lastB=",round(lastBM,1), ", bSm=",round(bSmM,2), 
#      ", lTFract=", round(locTrendFractM,2), sep='')				
  print(x$paramMeans)
  
#  cat(sprintf("NumTotalEvalEA: %d\n", x$numEvalEA))
#  cat(sprintf("NumTotalEvalLS: %d\n", x$numEvalLS))  
#  
#  ratio_effort <- x$numEvalEA/(x$numEvalEA+x$numEvalLS)
#  cat(sprintf("RatioEffort EA/LS: [%.0f/%.0f]\n", 100*ratio_effort, 100*(1-ratio_effort)))
#  
#  ratio_alg <- x$improvementEA/(x$improvementEA+x$improvementLS)
#  cat(sprintf("RatioImprovement EA/LS: [%.0f/%.0f]\n", 100*ratio_alg, 100*(1-ratio_alg)))
#  #cat(sprintf(("Restarts: %d\n", restarts))
#  
#  if((x$numTotalEA != 0) && (x$numTotalLS != 0)) {
#    cat(sprintf("PercentageNumImprovement[EA]: %d%%\n", round((x$numImprovementEA*100)/x$numTotalEA)))
#    cat(sprintf("PercentageNumImprovement[LS]: %d%%\n", round((x$numImprovementLS*100)/x$numTotalLS)))    
#  }
#  
#  cat(sprintf("Time[EA]: %.2f\n", x$timeMsEA))
#  cat(sprintf("Time[LS]: %.2f\n", x$timeMsLS))
#  cat(sprintf("Time[MA]: %.2f\n", x$timeMsMA))
#  cat(sprintf("RatioTime[EA/MA]: %.2f\n", 100*x$timeMsEA/x$timeMsMA))
#  cat(sprintf("RatioTime[LS/MA]: %.2f\n", 100*x$timeMsLS/x$timeMsMA))
#  
#  
#  cat("Fitness:\n",sep="")
#  print(x$fitness)
#  cat("Solution:\n",sep="") 
#  print(x$sol)
  
  invisible(x)
}












#  l=extract(samples)$l
#  lM=apply(l,2,mean)
#
#  b=extract(samples)$b
#  bM=apply(b,2,mean)
#  
#  levSm = extract(samples)$levSm
#  levSmM=mean(levSm)
#  
#  coefTrend = extract(samples)$coefTrend
#  coefTrendM=mean(coefTrend)
#  
#  powTrend = extract(samples)$powTrend
#  powTrendM=mean(powTrend)
#  
#  sigma = extract(samples)$sigma
#  sigmaM=mean(sigma)
#  
#  bSm = extract(samples)$bSm
#  bSmM=mean(bSm)
#  
#  locTrendFract = extract(samples)$locTrendFract
#  locTrendFractM=mean(locTrendFract)
#
#  if(isLGT) {
#    
#    nu = extract(samples)$nu
#    nuM=mean(nu)
#    
#    powx = extract(samples)$powx
#    powxM=mean(powx)
#    
#    offsetSigma=extract(samples)$offsetSigma
#    offsetSigmaM=mean(offsetSigma)
#  }

#  lastDetails=paste("coefT=",round(coefTrendM,2), ", powT=",round(powTrendM,2),
#      ", sigma=",round(sigmaM,2),", powx=",round(powxM,2),", offsetS=",round(offsetSigmaM,2), 
#      ", nu=",round(nuM,2), ", lSm=",round(levSmM,2), 
#      ", lastB=",round(lastBM,1), ", bSm=",round(bSmM,2), 
#      ", lTFract=", round(locTrendFractM,2), sep='')				
#  print(lastDetails)
#
#paste("lastLevel",round(lastLevelM,1))
#paste("levSm",levSmM)
#paste("nu",nuM)
#paste("coefTrend",coefTrendM)
#paste("powTrend",powTrendM)
#paste("sigma",sigmaM)
#paste("powx",powxM)
#paste("offsetSigma",round(offsetSigmaM,2))
#paste("lastB",round(lastBM,2))
#paste("bSm",round(bSmM,2))
#paste("locTrendFract",round(locTrendFractM,2))


#lastLevelM
#lastLevel
#powTrend
#coefTrend
#lastB
#levSm
#bSm
#sigma
#locTrendFract
#powx
#nu
#offsetSigma

#forc(lastLevelM = lastLevelM, NUM_OF_TRIALS = NUM_OF_TRIALS, h = h, mc_range = mc_range, 
#    lastLevel = lastLevel, powTrend = powTrend, coefTrend = coefTrend, lastB = lastB, 
#    levSm = levSm, bSm = bSm, sigma = sigma, locTrendFract = locTrendFract, isLGT = isLGT, 
#    powx = powx, nu = nu, offsetSigma = offsetSigma, MAX_VAL = MAX_VAL, MIN_VAL = MIN_VAL)

#  out <- list()
#  out[["yf"]] <- yf
#  out[["x"]] <- y.orig
#  class(out) <- "lgt"

#  out