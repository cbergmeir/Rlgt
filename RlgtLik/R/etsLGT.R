etsLGT <- function(y, model="ZZZ", 
    alpha=NULL, beta=NULL, gamma=NULL, phi=NULL, lambda=NULL, rho=NULL, additive.only=FALSE, 
    lower=c(rep(0.0001,3), 0, -1, 0), upper=c(rep(0.9999,3), 1, 1, 1),
    opt.crit=c("lik","amse","mse","sigma","mae"), nmse=3, bounds=c("both","usual","admissible"),
    ic=c("aicc","aic","bic"),restrict=TRUE, allow.multiplicative.trend=FALSE,
    use.initial.values=FALSE, solver="optim_c", ...)
{
  #lambdaBC=NULL,
  #biasadj=FALSE,
  #damped=NULL,
  
  damped <- TRUE
  model="AAN"

  #dataname <- substitute(y)
  opt.crit <- match.arg(opt.crit)
  bounds <- match.arg(bounds)
  ic <- match.arg(ic)

  if(any(class(y) %in% c("data.frame","list","matrix","mts")))
    stop("y should be a univariate time series")
  y <- as.ts(y)

#  # Check if data is constant
#  if (missing(model) & is.constant(y))
#    return(ses(y, alpha=0.99999, initial='simple')$model)

  # Remove missing values near ends
  ny <- length(y)
  y <- na.contiguous(y)
  if(ny != length(y))
    warning("Missing values encountered. Using longest contiguous portion of time series")

  orig.y <- y

  if(nmse < 1 | nmse > 30)
    stop("nmse out of range")
  m <- frequency(y)

  if(any(upper < lower))
    stop("Lower limits must be less than upper limits")

  errortype  <- substr(model,1,1)
  trendtype  <- substr(model,2,2)
  seasontype <- substr(model,3,3)

  if(!is.element(errortype,c("M","A","Z")))
    stop("Invalid error type")
  if(!is.element(trendtype,c("N","A","M","Z")))
    stop("Invalid trend type")
  if(!is.element(seasontype,c("N","A","M","Z")))
    stop("Invalid season type")

  if(m < 1 | length(y) <= m)
  {
    #warning("I can't handle data with frequency less than 1. Seasonality will be ignored.")
    seasontype <- "N"
  }
  if(m == 1)
  {
    if(seasontype=="A" | seasontype=="M")
      stop("Nonseasonal data")
    else
      substr(model,3,3) <- seasontype <- "N"
  }
  if(m > 24)
  {
    if(is.element(seasontype,c("A","M")))
      stop("Frequency too high")
    else if(seasontype=="Z")
    {
      warning("I can't handle data with frequency greater than 24. Seasonality will be ignored. Try stlf() if you need seasonal forecasts.")
      substr(model,3,3) <- seasontype <- "N"
      #m <- 1
    }
  }

  # Check inputs
  if(restrict)
  {
    if((errortype=="A" & (trendtype=="M" | seasontype=="M")) |
        (errortype=="M" & trendtype=="M" & seasontype=="A") |
        (additive.only & (errortype=="M" | trendtype=="M" | seasontype=="M")))
      stop("Forbidden model combination")
  }

  data.positive <- (min(y) > 0)

  if(!data.positive & errortype=="M")
    stop("Inappropriate model for data with negative or zero values")

  n <- length(y)

  # Fit model (assuming only one nonseasonal model)
  if(errortype=="Z")
    errortype <- c("A","M")
  if(trendtype=="Z")
  {
    if(allow.multiplicative.trend)
      trendtype <- c("N","A","M")
    else
      trendtype <- c("N","A")
  }
  if(seasontype=="Z")
    seasontype <- c("N","A","M")
  if(is.null(damped))
    damped <- c(TRUE,FALSE)
  best.ic <- Inf
  for(i in 1:length(errortype))
  {
    for(j in 1:length(trendtype))
    {
      for(k in 1:length(seasontype))
      {
        for(l in 1:length(damped))
        {
          if(trendtype[j]=="N" & damped[l])
            next
          if(restrict)
          {
            if(errortype[i]=="A" & (trendtype[j]=="M" | seasontype[k]=="M"))
              next
            if(errortype[i]=="M" & trendtype[j]=="M" & seasontype[k]=="A")
              next
            if(additive.only & (errortype[i]=="M" | trendtype[j]=="M" | seasontype[k]=="M"))
              next
          }
          if(!data.positive & errortype[i]=="M")
            next
          fit <- etsmodel(y,errortype[i],trendtype[j],seasontype[k],damped[l],alpha,beta,gamma,phi,lambda,rho,
              lower=lower,upper=upper,opt.crit=opt.crit,nmse=nmse,bounds=bounds, solver=solver, ...)
          fit.ic <- switch(ic,aic=fit$aic,bic=fit$bic,aicc=fit$aicc)
          if(!is.na(fit.ic))
          {
            if(fit.ic < best.ic)
            {
              model <- fit
              best.ic <- fit.ic
              best.e <- errortype[i]
              best.t <- trendtype[j]
              best.s <- seasontype[k]
              best.d <- damped[l]
            }
          }
        }
      }
    }
  }
  if(best.ic == Inf)
    stop("No model able to be fitted")

  model$m <- m
  model$method <- paste("ETS(",best.e,",",best.t,ifelse(best.d,"d",""),",",best.s,")",sep="")
  model$components <- c(best.e,best.t,best.s,best.d)
  model$call <- match.call()
  model$initstate <- model$states[1,]
  model$sigma2 <- mean(model$residuals^2,na.rm=TRUE)
  model$x <- orig.y

  return(structure(model,class="etsLGT"))
}




 myRequire <- function(libName) {

   req.suc <- require(libName, quietly=TRUE, character.only=TRUE)
   if(!req.suc) stop("The ",libName," package is not available.")

   req.suc
 }

 getNewBounds <- function(par, lower, upper, nstate) {

   myLower <- NULL
   myUpper <- NULL

   if("alpha" %in% names(par)) {
     myLower <- c(myLower, lower[1])
     myUpper <- c(myUpper, upper[1])
   }
   if("beta" %in% names(par)) {
     myLower <- c(myLower, lower[2])
     myUpper <- c(myUpper, upper[2])
   }
   if("gamma" %in% names(par)) {
     myLower <- c(myLower, lower[3])
     myUpper <- c(myUpper, upper[3])
   }
   if("phi" %in% names(par)) {
     myLower <- c(myLower, lower[4])
     myUpper <- c(myUpper, upper[4])
   }
   if("lambda" %in% names(par)) {
     myLower <- c(myLower, lower[5])
     myUpper <- c(myUpper, upper[5])
   }
   if("rho" %in% names(par)) {
     myLower <- c(myLower, lower[6])
     myUpper <- c(myUpper, upper[6])
   }
   
   myLower <- c(myLower,rep(-1e8,nstate))
   myUpper <- c(myUpper,rep(1e8,nstate))

   list(lower=myLower, upper=myUpper)
 }


etsmodel <- function(y, errortype, trendtype, seasontype, damped,
    alpha=NULL, beta=NULL, gamma=NULL, phi=NULL, lambda=NULL, rho=NULL,
    lower, upper, opt.crit, nmse, bounds, maxit=2000, control=NULL, seed=NULL, trace=FALSE, solver="optim_c")
{

  tsp.y <- tsp(y)
  if(is.null(tsp.y))
    tsp.y <- c(1,length(y),1)
  if(seasontype != "N")
    m <- tsp.y[3]
  else
    m <- 1

  # Initialize smoothing parameters
  par <- initparam(alpha,beta,gamma,phi,lambda, rho, trendtype,seasontype,lower,upper,m)
  names(alpha) <- names(beta) <- names(gamma) <- names(phi) <- names(lambda) <- names(rho) <- NULL
  par.noopt <- c(alpha=alpha,beta=beta,gamma=gamma,phi=phi, lambda=lambda, rho=rho)
  if(!is.null(par.noopt))
    par.noopt <- c(na.omit(par.noopt))
  if(!is.na(par["alpha"]))
    alpha <- par["alpha"]
  if(!is.na(par["beta"]))
    beta <- par["beta"]
  if(!is.na(par["gamma"]))
    gamma <- par["gamma"]
  if(!is.na(par["phi"]))
    phi <- par["phi"]
  if(!is.na(par["lambda"]))
    lambda <- par["lambda"]
  if(!is.na(par["rho"]))
    rho <- par["rho"]
  
#    if(errortype=="M" | trendtype=="M" | seasontype=="M")
#        bounds="usual"
  if(!check.param(alpha,beta,gamma,phi,lambda,rho,lower,upper,bounds,m))
  {
    #print(paste("Model: ETS(",errortype,",",trendtype,ifelse(damped,"d",""),",",seasontype,")",sep=""))
    stop("Parameters out of range")
  }

  # Initialize state
  init.state <- initstate(y,trendtype,seasontype)
  nstate <- length(init.state)
  par <- c(par,init.state)
  lower <- c(lower,rep(-Inf,nstate))
  upper <- c(upper,rep(Inf,nstate))

  np <- length(par)
  if(np >= length(y)-1) # Not enough data to continue
    return(list(aic=Inf,bic=Inf,aicc=Inf,mse=Inf,amse=Inf,fit=NULL,par=par,states=init.state))

#-------------------------------------------------

  if(is.null(seed)) seed <- 1000*runif(1)

   if(solver=="malschains_c") {

     malschains <- NULL
     if(!myRequire("Rmalschains"))
       stop("malschains optimizer unavailable")

     func <- NULL
     env <- NULL

       env <- etsTargetFunctionInit(par=par, y=y, nstate=nstate, errortype=errortype, trendtype=trendtype,
           seasontype=seasontype, damped=damped, par.noopt=par.noopt, lowerb=lower, upperb=upper,
           opt.crit=opt.crit, nmse=nmse, bounds=bounds, m=m,pnames=names(par),pnames2=names(par.noopt))

       func <- .Call("etsGetTargetFunctionRmalschainsPtr", PACKAGE="RlgtLik")

     myBounds <- getNewBounds(par, lower, upper, nstate)

     if(is.null(control)) {
       control <- Rmalschains::malschains.control(ls="simplex", lsOnly=TRUE)
     }

     control$optimum <- if(opt.crit=="lik") -1e12 else 0

     fredTmp <- Rmalschains::malschains(func, env=env, lower=myBounds$lower, upper=myBounds$upper,
         maxEvals=maxit, seed=seed, initialpop=par, control=control)

     fred <- NULL
     fred$par <- fredTmp$sol

     fit.par <- fred$par

     names(fit.par) <- names(par)

  } else if(solver=="optim_c"){

    env <- etsTargetFunctionInit(par=par, y=y, nstate=nstate, errortype=errortype, trendtype=trendtype,
        seasontype=seasontype, damped=damped, par.noopt=par.noopt, lowerb=lower, upperb=upper,
        opt.crit=opt.crit, nmse=as.integer(nmse), bounds=bounds, m=m,pnames=names(par),pnames2=names(par.noopt))

    fred <- .Call("etsNelderMead", par, env, -Inf,
        sqrt(.Machine$double.eps), 1.0, 0.5, 2.0, trace, maxit, PACKAGE="RlgtLik")

    fit.par <- fred$par

    names(fit.par) <- names(par)

   }


#-------------------------------------------------

  init.state <- fit.par[(np-nstate+1):np]
  # Add extra state
  if(seasontype!="N")
    init.state <- c(init.state, m*(seasontype=="M") - sum(init.state[(2+(trendtype!="N")):nstate]))

  if(!is.na(fit.par["alpha"]))
    alpha <- fit.par["alpha"]
  if(!is.na(fit.par["beta"]))
    beta <- fit.par["beta"]
  if(!is.na(fit.par["gamma"]))
    gamma <- fit.par["gamma"]
  if(!is.na(fit.par["phi"]))
    phi <- fit.par["phi"]
  if(!is.na(fit.par["lambda"]))
    lambda <- fit.par["lambda"]
  if(!is.na(fit.par["rho"]))
    rho <- fit.par["rho"]
  
  e <- pegelsresid.C(y,m,init.state,errortype,trendtype,seasontype,alpha,beta,gamma,phi,lambda,rho,nmse)

  np <- np + 1
  ny <- length(y)
  aic <- e$lik + 2*np
  bic <- e$lik + log(ny)*np
  aicc <- aic +  2*np*(np+1)/(ny-np-1)

  mse <- e$amse[1]
  amse <- mean(e$amse)

  states=ts(e$states,frequency=tsp.y[3],start=tsp.y[1]-1/tsp.y[3])
  colnames(states)[1] <- "l"
  if(trendtype!="N")
    colnames(states)[2] <- "b"
  if(seasontype!="N")
    colnames(states)[(2+(trendtype!="N")):ncol(states)] <- paste("s",1:m,sep="")

  #tmp <- c("alpha",rep("beta",trendtype!="N"),rep("gamma",seasontype!="N"),rep("phi",damped))
  fit.par <- c(fit.par,par.noopt)
#    fit.par <- fit.par[order(names(fit.par))]
  if(errortype=="A")
    fits <- y-e$e
  else
    fits <- y/(1+e$e)

  return(list(loglik=-0.5*e$lik,aic=aic,bic=bic,aicc=aicc,mse=mse,amse=amse,fit=fred,residuals=ts(e$e,frequency=tsp.y[3],start=tsp.y[1]),fitted=ts(fits,frequency=tsp.y[3],start=tsp.y[1]),
          states=states,par=fit.par))
}


etsTargetFunctionInit <- function(par,y,nstate,errortype,trendtype,seasontype,damped,par.noopt,lowerb,upperb,
    opt.crit,nmse,bounds,m,pnames,pnames2)
{

  names(par) <- pnames
  names(par.noopt) <- pnames2
  alpha <- c(par["alpha"],par.noopt["alpha"])["alpha"]
  if(is.na(alpha))
    stop("alpha problem!")
  if(trendtype!="N")
  {
    beta <- c(par["beta"],par.noopt["beta"])["beta"]
    lambda <- c(par["lambda"],par.noopt["lambda"])["lambda"]
    rho <- c(par["rho"],par.noopt["rho"])["rho"]
    if(is.na(beta))
      stop("beta Problem!")
  }
  else
    beta <- NULL
  if(seasontype!="N")
  {
    gamma <- c(par["gamma"],par.noopt["gamma"])["gamma"]
    if(is.na(gamma))
      stop("gamma Problem!")
  }
  else
  {
    m <- 1
    gamma <- NULL
  }

  phi <- c(par["phi"],par.noopt["phi"])["phi"]
    if(is.na(phi))
      stop("phi Problem!")

  #determine which values to optimize and which ones are given by the user/not needed
  optAlpha <- !is.null(alpha)
  optBeta <- !is.null(beta)
  optGamma <- !is.null(gamma)
  optPhi <- !is.null(phi)
  optLambda <- !is.null(lambda)
  optRho <- !is.null(rho)
  
  givenAlpha <- FALSE
  givenBeta <- FALSE
  givenGamma <- FALSE
  givenPhi <- FALSE
  givenLambda <- FALSE
  givenRho <- FALSE
  
  if(!is.null(par.noopt["alpha"])) if(!is.na(par.noopt["alpha"])) {
      optAlpha <- FALSE
      givenAlpha <- TRUE
    }
  if(!is.null(par.noopt["beta"])) if(!is.na(par.noopt["beta"])) {
      optBeta <- FALSE
      givenBeta <- TRUE
    }
  if(!is.null(par.noopt["gamma"])) if(!is.na(par.noopt["gamma"])) {
      optGamma <- FALSE
      givenGamma <- TRUE
    }
  if(!is.null(par.noopt["phi"])) if(!is.na(par.noopt["phi"])) {
      optPhi <- FALSE
      givenPhi <- TRUE
    }
  if(!is.null(par.noopt["lambda"])) if(!is.na(par.noopt["lambda"])) {
      optLambda <- FALSE
      givenLambda <- TRUE
    }
  if(!is.null(par.noopt["rho"])) if(!is.na(par.noopt["rho"])) {
      optRho <- FALSE
      givenRho <- TRUE
    }
  
  if(trendtype == "N")
    beta <- 0;
  if(seasontype == "N")
    gamma <- 0;

  #lambda <- phi
  #rho <- 0
  
#  cat("alpha: ", alpha)
#  cat(" beta: ", beta)
#  cat(" gamma: ", gamma)
#  cat(" phi: ", phi, "\n")
#
#  cat("useAlpha: ", useAlpha)
#  cat(" useBeta: ", useBeta)
#  cat(" useGamma: ", useGamma)
#  cat(" usePhi: ", usePhi, "\n")

  env <- new.env()

  res <- .Call("etsTargetFunctionInit", y=y, nstate=nstate, errortype=switch(errortype,"A"=1,"M"=2),
      trendtype=switch(trendtype,"N"=0,"A"=1,"M"=2), seasontype=switch(seasontype,"N"=0,"A"=1,"M"=2),
      damped=damped, lowerb=lowerb, upperb=upperb,
      opt.crit=opt.crit, nmse=as.integer(nmse), bounds=bounds, m=m,
      optAlpha, optBeta, optGamma, optPhi, optLambda, optRho,
      givenAlpha, givenBeta, givenGamma, givenPhi, givenLambda, givenRho,
      alpha, beta, gamma, phi, lambda, rho, env, PACKAGE="RlgtLik")
  res
}



initparam <- function(alpha,beta,gamma,phi,lambda,rho,trendtype,seasontype,lower,upper,m)
{
  if(any(lower > upper))
    stop("Inconsistent parameter boundaries")

  # Select alpha
  if(is.null(alpha))
  {
    alpha <- lower[1] + 0.5*(upper[1]-lower[1])/m
    par <- c(alpha=alpha)
  }
  else
    par <- numeric(0)

  # Select beta
  if(trendtype !="N" & is.null(beta))
  {
    # Ensure beta < alpha
    upper[2] <- min(upper[2], alpha)
    beta <- lower[2] + 0.1*(upper[2]-lower[2])
    par <- c(par,beta=beta)
  }

  # Select gamma
  if(seasontype != "N" & is.null(gamma))
  {
    # Ensure gamma < 1-alpha
    upper[3] <- min(upper[3], 1-alpha)
    gamma <- lower[3] + 0.05*(upper[3]-lower[3])
    par <- c(par,gamma=gamma)
  }

  # Select phi
  if(is.null(phi))
  {
    phi <- lower[4] + .99*(upper[4]-lower[4])
    par <- c(par,phi=phi)
  }

  if(is.null(lambda))
  {
    lambda <- lower[5] + .99*(upper[5]-lower[5])
    par <- c(par,lambda=lambda)
  }

  if(is.null(rho))
  {
    rho <- 0 #lower[4] + .99*(upper[4]-lower[4])
    par <- c(par,rho=rho)
  }
  
  return(par)
}

check.param <- function(alpha,beta,gamma,phi,lambda,rho,lower,upper,bounds,m)
{
#  if(bounds != "admissible")
#  {
    if(!is.null(alpha))
    {
      if(alpha < lower[1] | alpha > upper[1])
        return(0)
    }
    if(!is.null(beta))
    {
      if(beta < lower[2] | beta > alpha | beta > upper[2])
        return(0)
    }
    if(!is.null(phi))
    {
      if(phi < lower[4] | phi > upper[4])
        return(0)
    }
    if(!is.null(lambda))
    {
      if(lambda < lower[5] | lambda > upper[5])
        return(0)
    }
    if(!is.null(rho))
    {
      if(rho < lower[6] | rho > upper[6])
        return(0)
    }
    if(!is.null(gamma))
    {
      if(gamma < lower[3] | gamma > 1-alpha | gamma > upper[3])
        return(0)
    }
#  }
#  if(bounds != "usual")
#  {
#    if(!admissible(alpha,beta,gamma,phi,m))
#      return(0)
#  }
  return(1)
}

initstate <- function(y,trendtype,seasontype)
{
#  if(seasontype!="N")
#  {
#    # Do decomposition
#    m <- frequency(y)
#    n <- length(y)
#    if(n < 4)
#      stop("You've got to be joking (not enough data).")
#    else if(n < 3*m) # Fit simple Fourier model.
#    {
#      fouriery <- fourier(y,1)
#      fit <- tslm(y ~ trend + fouriery)
#      if(seasontype=="A")
#        y.d <- list(seasonal=y -fit$coef[1] - fit$coef[2]*(1:n))
#      else # seasontype=="M". Biased method, but we only need a starting point
#        y.d <- list(seasonal=y / (fit$coef[1] + fit$coef[2]*(1:n)))
#    }
#    else # n is large enough to do a decomposition
#      y.d <- decompose(y,type=switch(seasontype, A="additive", M="multiplicative"))
#
#    init.seas <- rev(y.d$seasonal[2:m]) # initial seasonal component
#    names(init.seas) <- paste("s",0:(m-2),sep="")
#    # Seasonally adjusted data
#    if(seasontype=="A")
#      y.sa <- y-y.d$seasonal
#    else
#    {
#      init.seas <- pmax(init.seas, 1e-2) # We do not want negative seasonal indexes
#      if(sum(init.seas) > m)
#      	init.seas <- init.seas/sum(init.seas + 1e-2)
#      y.sa <- y/pmax(y.d$seasonal, 1e-2)
#    }
#  }
#  else # non-seasonal model
#  {
    m <- 1
    init.seas <- NULL
    y.sa <- y
#  }

  maxn <- min(max(10,2*m),length(y.sa))

#  if(trendtype=="N")
#  {
#    l0 <- mean(y.sa[1:maxn])
#    b0 <- NULL
#  }
#  else  # Simple linear regression on seasonally adjusted data
#  {
    fit <- lsfit(1:maxn,y.sa[1:maxn])
#    if(trendtype=="A")
#    {
      l0 <- fit$coef[1]
      b0 <- fit$coef[2]
      # If error type is "M", then we don't want l0+b0=0.
      # So perturb just in case.
      if(abs(l0+b0) < 1e-8)
      {
        l0 <- l0*(1+1e-3)
        b0 <- b0*(1-1e-3)
      }
#    }
#    else #if(trendtype=="M")
#    {
#      l0 <- fit$coef[1]+fit$coef[2] # First fitted value
#      if(abs(l0) < 1e-8)
#        l0 <- 1e-7
#      b0 <- (fit$coef[1] + 2*fit$coef[2])/l0 # Ratio of first two fitted values
#      l0 <- l0/b0 # First fitted value divided by b0
#      if(abs(b0) > 1e10) # Avoid infinite slopes
#        b0 <- sign(b0)*1e10
#      if(l0 < 1e-8 | b0 < 1e-8) # Simple linear approximation didn't work.
#      {
#        l0 <- max(y.sa[1],1e-3)
#        b0 <- max(y.sa[2]/y.sa[1],1e-3)
#      }
#    }
#  }

  names(l0) <- "l"
  if(!is.null(b0))
    names(b0) <- "b"
  return(c(l0,b0,init.seas))
}


pegelsresid.C <- function(y,m,init.state,errortype,trendtype,seasontype,alpha,beta,gamma,phi,lambda,rho,nmse)
{
  n <- length(y)
  p <- length(init.state)
  x <- numeric(p*(n+1))
  x[1:p] <- init.state
  e <- numeric(n)
  lik <- 0;
#  if(!damped)
#    phi <- 1;
  if(trendtype == "N")
    beta <- 0;
  if(seasontype == "N")
    gamma <- 0;

  amse <- numeric(nmse)

  Cout <- .C("etscalc",
      as.double(y),
      as.integer(n),
      as.double(x),
      as.integer(m),
      as.integer(switch(errortype,"A"=1,"M"=2)),
      as.integer(switch(trendtype,"N"=0,"A"=1,"M"=2)),
      as.integer(switch(seasontype,"N"=0,"A"=1,"M"=2)),
      as.double(alpha),
      as.double(beta),
      as.double(gamma),
      as.double(phi),
      as.double(lambda),
      as.double(rho),
      as.double(e),
      as.double(lik),
      as.double(amse),
      as.integer(nmse),
      PACKAGE="RlgtLik")
  if(!is.na(Cout[[15]]))
  {
    if(abs(Cout[[15]]+99999) < 1e-7)
      Cout[[15]] <- NA
  }
  tsp.y <- tsp(y)
  e <- ts(Cout[[14]])
  tsp(e) <- tsp.y

  return(list(lik=Cout[[15]], amse=Cout[[16]], e=e, states=matrix(Cout[[3]], nrow=n+1, ncol=p, byrow=TRUE)))
}

#admissible <- function(alpha,beta,gamma,phi,lambda,rho,m)
#{
#  if(is.null(phi))
#    phi <- 1
#  if(phi < 0 | phi > 1+1e-8)
#    return(0)
#  if(is.null(gamma))
#  {
#    if(alpha < 1-1/phi | alpha > 1+1/phi)
#      return(0)
#    if(!is.null(beta))
#    {
#      if(beta < alpha * (phi-1) | beta > (1+phi)*(2-alpha))
#        return(0)
#    }
#  }
#  else if(m > 1) # Seasonal model
#  {
#    if(is.null(beta))
#      beta <- 0
#    if(gamma < max(1-1/phi-alpha,0) | gamma > 1+1/phi-alpha)
#      return(0)
#    if(alpha < 1-1/phi-gamma*(1-m+phi+phi*m)/(2*phi*m))
#      return(0)
#    if(beta < -(1-phi)*(gamma/m+alpha))
#      return(0)
#
#    # End of easy tests. Now use characteristic equation
#    P <- c(phi*(1-alpha-gamma),alpha+beta-alpha*phi+gamma-1,rep(alpha+beta-alpha*phi,m-2),(alpha+beta-phi),1)
#    roots <- polyroot(P)
#
#    #cat("maxpolyroots: ", max(abs(roots)), "\n")
#
#    if(max(abs(roots)) > 1+1e-10)
#      return(0)
#  }
#  # Passed all tests
#  return(1)
#}





print.etsLGT <- function(x,...)
{
  cat(paste(x$method, "\n\n"))
  cat(paste("Call:\n", deparse(x$call), "\n\n"))
  ncoef <- length(x$initstate)
  if(!is.null(x$lambda))
    cat("  Box-Cox transformation: lambda=",round(x$lambda,4), "\n\n")
  
  cat("  Smoothing parameters:\n")
  cat(paste("    alpha =", round(x$par["alpha"], 4), "\n"))
  if(x$components[2]!="N")
    cat(paste("    beta  =", round(x$par["beta"], 4), "\n"))
  if(x$components[3]!="N")
    cat(paste("    gamma =", round(x$par["gamma"], 4), "\n"))
  if(x$components[4]!="FALSE")
    cat(paste("    phi   =", round(x$par["phi"], 4), "\n"))
  
  cat("\n  Initial states:\n")
  cat(paste("    l =", round(x$initstate[1], 4), "\n"))
  if (x$components[2]!="N")
    cat(paste("    b =", round(x$initstate[2], 4), "\n"))
  else
  {
    x$initstate <- c(x$initstate[1], NA, x$initstate[2:ncoef])
    ncoef <- ncoef+1
  }
  if (x$components[3]!="N")
  {
    cat("    s=")
    if (ncoef <= 8)
      cat(round(x$initstate[3:ncoef], 4))
    else
    {
      cat(round(x$initstate[3:8], 4))
      cat("\n           ")
      cat(round(x$initstate[9:ncoef], 4))
    }
    cat("\n")
  }
  
  cat("\n  sigma:  ")
  cat(round(sqrt(x$sigma2),4))
  if(!is.null(x$aic))
  {
    stats <- c(x$aic,x$aicc,x$bic)
    names(stats) <- c("AIC","AICc","BIC")
    cat("\n\n")
    print(stats)
  }
#    cat("\n  AIC:    ")
#    cat(round(x$aic,4))
#    cat("\n  AICc:   ")
#    cat(round(x$aicc,4))
#    cat("\n  BIC:    ")
#    cat(round(x$bic,4))
}


### PLOT COMPONENTS
plot.etsLGT <- function(x,...)
{
  if(!is.null(x$lambda))
    y <- BoxCox(x$x,x$lambda)
  else
    y <- x$x
  if(x$components[3]=="N" & x$components[2]=="N")
  {
    plot(cbind(observed=y, level=x$states[,1]),
        main=paste("Decomposition by",x$method,"method"),...)
  }
  else if(x$components[3]=="N")
  {
    plot(cbind(observed=y, level=x$states[,1], slope=x$states[,"b"]),
        main=paste("Decomposition by",x$method,"method"),...)
  }
  else if(x$components[2]=="N")
  {
    plot(cbind(observed=y, level=x$states[,1], season=x$states[,"s1"]),
        main=paste("Decomposition by",x$method,"method"),...)
  }
  else
  {
    plot(cbind(observed=y, level=x$states[,1], slope=x$states[,"b"],
            season=x$states[,"s1"]),
        main=paste("Decomposition by",x$method,"method"),...)
  }
}

summary.etsLGT <- function(object,...)
{
  print(object)
  cat("\nTraining set error measures:\n")
  print(accuracy(object))
}

coef.etsLGT <- function(object,...)
{
  object$par
}

fitted.etsLGT <- function(object, h=1, ...){
  if(h==1){
    return(object$fitted)
  }
  else{
    return(hfitted(object=object, h=h, FUN="ets", ...))
  }
}

logLik.etsLGT <- function(object,...)
{
  structure(object$loglik,df=length(object$par),class="logLik")
}

is.etsLGT <- function(x){
  inherits(x, "etsLGT")
}
